function ihMT_R1vsM0b_correlation_3prot(obj)
%% Generate the M0B mapping to R1 from simulation results and acquired data

%  NOTE: niak_read_vol replaced with minc_read  %
%  NOTE: niak_write_vol replaced with minc_write %


% Where outputs will be saved/ previous results are 
OutputDir = obj.options.R1vsM0bMapping_OutputDirectory;
%OutputDir = '/Users/amiedemmans/Desktop/GitHub/qMRLab_test_dir/Test1';
%% Load images:

% Where Images are you will be using 
DATADIR = obj.options.R1vsM0bMapping_DataDirectory;
%DATADIR = '/Users/amiedemmans/Desktop/GitHub/Image_files/ihMT_images_from_chris/';

%image names:
mtw_fn = {'dual_reg.mnc' 'pos_reg.mnc' 'neg_reg.mnc' 'sparseMP2RAGE_M0.mnc' 'sparseMP2RAGE_T1.mnc'};

for i = 1:size(mtw_fn,2)
    fn = fullfile(DATADIR,mtw_fn{i});
    [hdr, img] = minc_read(fn); % niak_read_vol
    comb_mtw(:,:,:,i) = img; %.img;
end


%% Load the mask
fn = fullfile(DATADIR,'itkmask.mnc');
[~, mask] = minc_read(fn);
mask1 = permute(mask,[2 3 1]); % conversion between minc and nii reorients it

% If making a mask 
% mask = zeros(size(lfa)); 
% mask(lfa>175) = 1;

%% Some B1 issues so lets try and load that
[~, b1] = minc_read(fullfile(DATADIR,'resampled_b1field.mnc')); 
b1 = double(b1);
b1 = limitHandler(b1, 0.5,1.4);
b1 = ihMT_imgaussfilt3_withMask(b1, mask1, 5); %light smoothing to the map

figure; imshow3Dfull(b1, [0.6 1.2],jet)

%% List the images 
dual = comb_mtw(:,:,:,1);
pos = comb_mtw(:,:,:,2);
neg = comb_mtw(:,:,:,3);

spApp_mp2 = comb_mtw(:,:,:,4); 
spT1_map = comb_mtw(:,:,:,5); 


%% Mask -> bet result touched up in itk, then threshold CSF and some dura
% Does this go away because mp2rage section is gone 
%mask = mask1;
mask(spT1_map > 2500) = 0;
mask(spT1_map < 650) = 0;
mask(isnan(spT1_map)) = 0;
mask = bwareaopen(mask, 10000,6);
figure; imshow3Dfullseg(spT1_map, [300 2500],mask)


%% Compute MTsat ihMTsat

%% Protocol 
flipA = obj.Prot.PulseSequenceParams.Mat(3);
TR = obj.Prot.PulseSequenceParams.Mat(4);
DummyEcho = obj.Prot.PulseSequenceParams.Mat(10);
echoSpacing = obj.Prot.PulseSequenceParams.Mat(14);
numExcitation = obj.Prot.PulseSequenceParams.Mat(6) + DummyEcho;

sat_dual1 = ihMT_calcMTsatThruLookupTablewithDummyV3( dual, b1, spT1_map, mask1,spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_pos1  = ihMT_calcMTsatThruLookupTablewithDummyV3( pos, b1, spT1_map, mask1, spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_neg1  = ihMT_calcMTsatThruLookupTablewithDummyV3( neg, b1, spT1_map, mask1, spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% figure; imshow3Dfull(sat_dual1 , [0 0.06], jet); figure; imshow3Dfull(sat_pos1 , [0 0.05], jet);  figure; imshow3Dfull(sat_neg1 , [0 0.05], jet); 

%% With MTsat maps made, perform M0b mapping

% load in the fit results for VFA - Optimal
fitValues_S = load(fullfile(OutputDir,'fitValues_S.mat'));
fitValues_S = fitValues_S.fitValues;
fitValues_D = load(fullfile(OutputDir,'fitValues_D.mat'));
fitValues_D = fitValues_D.fitValues;


% need to convert to 1/s from 1/ms -> ONLY USE MP2RAGE values, VFA are too
% far off.
R1_s = (1./spT1_map) *1000;

% initialize matrices
M0b_dual = zeros(size(sat_dual1));
M0b_pos = zeros(size(sat_dual1));
M0b_neg = zeros(size(sat_dual1));


%% SPEED IT UP BY DOING A FEW AXIAL SLICES
axialStart = 126; % 65
axialStop = axialStart+3;%115;
% check
figure; imshow3Dfull(sat_dual1(:,axialStart:axialStop,:) , [0 0.06], jet)

b1_1 = obj.Prot.PulseSequenceParams.Mat(8);

tic %  
for i = 1:size(sat_dual1,1) % went to 149
    
    for j = axialStart:axialStop % 1:size(sat_dual,2) % for axial slices
        for k =  1:size(sat_dual1,3) % sagital slices  65
            
            if mask(i,j,k) > 0 %&& dual_s(i,j,k,3) > 0
                                
                 [M0b_dual(i,j,k), ~,  ~]  = ihMT_fit_M0b_v2( b1_1*b1(i,j,k), R1_s(i,j,k), sat_dual1(i,j,k), fitValues_D);
                 [M0b_pos(i,j,k),  ~,  ~]  = ihMT_fit_M0b_v2( b1_1*b1(i,j,k), R1_s(i,j,k), sat_pos1(i,j,k), fitValues_S);               
                 [M0b_neg(i,j,k),  ~,  ~]  = ihMT_fit_M0b_v2( b1_1*b1(i,j,k), R1_s(i,j,k), sat_neg1(i,j,k), fitValues_S);
                 
            end
        end
    end
    disp(i)
end
toc %% this took 30hours for 1mm isotropic full brain dataset. * was running fitting in another matlab
    % instance, so could be easily sped up running on its own and/or adding
    % the parfor loop. 


figure; imshow3Dfull(M0b_pos, [0 0.15],jet)
    
figure; imshow3Dfull(M0b_dual, [0 0.15],jet)

% export
%mkdir(fullfile(OutputDir,'processing'))
%hdr.file_name = fullfile(OutputDir,'processing/M0b_dual.mnc.gz'); 
minc_write('M0b_dual.mnc.gz', hdr, M0b_dual); % Need to add filename for each minc_write %
%hdr.file_name = fullfile(OutputDir,'processing/M0b_pos.mnc.gz');
minc_write('M0b_pos.mnc.gz', hdr, M0b_pos);
%hdr.file_name = fullfile(OutputDir,'processing/M0b_neg.mnc.gz'); 
minc_write('M0b_neg.mnc.gz', hdr, M0b_neg);


%% With M0B maps made, correlate with R1 and update the fitValues file. 

% use this fake mask to get rid of dura. 
tempMask = mask;
tempMask = imerode(tempMask, strel('sphere',2));
figure; imshow3Dfullseg(M0b_200_dual, [0 0.15],tempMask)

mkdir(fullfile(OutputDir,'figures'));

% Optimized Approach
fitValues_D  = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_8_dual, tempMask, fitValues_D, fullfile(OutputDir,'figures/R1vsM0b_8_dual.png'), fullfile(OutputDir,'fitValues_D_8.mat'));
fitValues_SP = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_8_pos, tempMask, fitValues_S, fullfile(OutputDir,'figures/R1vsM0b_8_pos.png'), fullfile(OutputDir,'fitValues_SP_8.mat'));
fitValues_SN = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_8_neg, tempMask, fitValues_S, fullfile(OutputDir,'figures/R1vsM0b_8_neg.png'), fullfile(OutputDir,'fitValues_SN_8.mat'));


%% Now use these results to B1 correct the data:
OutputDir = DATADIR;

b1_1 = obj.Prot.PulseSequenceParams.Mat(8);

corr_prot1_d = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_D);
corr_prot1_p = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_SP);
corr_prot1_n = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_SN);


% Part 2, apply correction map
sat_dual1_c = (sat_dual1 + sat_dual1.* corr_prot1_d) .* mask1;
sat_pos1_c  = (sat_pos1 + sat_pos1.* corr_prot1_p) .* mask1;
sat_neg1_c  = (sat_neg1 + sat_neg1.* corr_prot1_n) .* mask1;
ihmt1_c      = sat_dual1_c - (sat_pos1_c + sat_neg1_c)/2;


ihmt1_c = double(limitHandler(ihmt1_c,0, 0.05));

ihmt1_c( ihmt1_c >= 0.05) = 0;

%% View results
figure; imshow3Dfull(sat_dual1_c , [0 0.06], jet); 
figure; imshow3Dfull(sat_pos1_c , [0 0.06], jet)
figure; imshow3Dfull(sat_neg1_c , [0 0.06], jet); 

figure; imshow3Dfull(ihmt1_c , [0 0.06], jet)

%% Other things, save if you want

%hdr.file_name = strcat(DATADIR,'matlab/MTsat_dual_1.mnc.gz'); 
minc_write('MTsat_dual_1.mnc.gz', hdr, sat_dual1_c);
%hdr.file_name = strcat(DATADIR,'matlab/MTsat_pos_1.mnc.gz');
minc_write('MTsat_pos_1.mnc.gz', hdr, sat_pos1_c);
%hdr.file_name = strcat(DATADIR,'matlab/MTsat_neg_1.mnc.gz'); 
minc_write('MTsat_neg_1.mnc.gz', hdr, sat_neg1_c);
%hdr.file_name = strcat(DATADIR,'matlab/ihMTsat_1.mnc.gz');
minc_write('ihMTsat_1.mnc.gz', hdr, ihmt1_c);

%hdr.file_name = strcat(DATADIR,'matlab/b1.mnc.gz'); 
minc_write('b1.mnc.gz', hdr, b1);



%% 

ihmtSlice1 = squeeze( ihmt1_c(:,126,:));

figure; imagesc(ihmtSlice1); axis image;
colormap(gray)
caxis([0 0.03])
hold on
line([0,175], [165,166], 'Color', 'r');


ihmtProf1 = ihmtSlice1(165,:);

% normalize:
ihmtProf1 = ihmtProf1/ max(ihmtProf1);

figure
plot(ihmtProf1,'LineWidth',2); 
title('Line Profile (L-R)');
xlim([15, 160])
hold off
xlabel('Voxel Index');
ylabel('Relative ihMT_{sat}')
ax = gca; ax.FontSize = 20; 
set(gcf,'position',[10,400,1000,600])   



%% Stuff that's not needed 
%% Run MP-PCA denoising
% comb_mtw = double(comb_mtw);
% all_PCAcorr = MPdenoising(comb_mtw); % MPdenoising already exists 


%% separate the images then average the MTw ones

% dual = all_PCAcorr(:,:,:,1);
% pos  = all_PCAcorr(:,:,:,2);
% neg = all_PCAcorr(:,:,:,3);
% noMT = all_PCAcorr(:,:,:,4);
% 
% dual2 = all_PCAcorr(:,:,:,5);
% pos2  = all_PCAcorr(:,:,:,6);
% neg2 = all_PCAcorr(:,:,:,7);
% noMT2 = all_PCAcorr(:,:,:,8);
% 
% dual3 = all_PCAcorr(:,:,:,9);
% pos3  = all_PCAcorr(:,:,:,10);
% neg3 = all_PCAcorr(:,:,:,11);
% noMT3 = all_PCAcorr(:,:,:,12);
% 
% sp_mp2r_inv1  = all_PCAcorr(:,:,:,13);
% sp_mp2r_inv2 = all_PCAcorr(:,:,:,14);
% sp_mp2r_uni  = all_PCAcorr(:,:,:,15);


%% Now from the MP2RAGE:
  
% MP2RAGE.B0 = 3;           % in Tesla
% MP2RAGE.TR = 5;           % MP2RAGE TR in seconds
% MP2RAGE.TRFLASH = 6.4e-3; % TR of the GRE readout
% MP2RAGE.TIs = [0.94 2.83];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
% MP2RAGE.NZslices = [ceil(175/2) floor(175/2)];%  should be two values, number of excitations before k-space center, and number after. [Slices Per Slab * [PartialFourierInSlice-0.5  0.5] ]
% MP2RAGE.FlipDegrees = [4 5];% Flip angle of the two readouts in degrees
% 
% 
% MP2RAGEimg.img = sp_mp2r_uni; % load_untouch_nii(MP2RAGE.filenameUNI);
% MP2RAGEINV2img.img = sp_mp2r_inv2; % load_untouch_nii(MP2RAGE.filenameINV2);
% %B1.img = permute(b1,[3 2 1]);
% B1.img = b1;
% brain.img = mask;
% 
% tic
% [ spT1map, spMP2RAGEcorrected, spAppmap2] = ihMT_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);
% toc
% 
% spT1_map = spT1map.img;
% spApp_mp2 = spAppmap2.img ;
% 
% spT1_map = limitHandler(spT1_map);
% spApp_mp2 = double(limitHandler(spApp_mp2));
% 
% figure; imshow3Dfull(spT1_map, [300 2500],jet)
% figure; imshow3Dfull(spApp_mp2 , [00 6000])
% 
% 
% minc_write('sparseMP2RAGE_T1.mnc', hdr, spT1_map);
% minc_write('sparseMP2RAGE_M0.mnc', hdr, spApp_mp2);
% hdr.file_name = fullfile(DATADIR,'matlab/sparseMP2RAGE_T1.mnc.gz'); niak_write_vol(hdr, spT1_map);
% hdr.file_name = fullfile(DATADIR,'matlab/sparseMP2RAGE_M0.mnc.gz'); niak_write_vol(hdr, spApp_mp2);










