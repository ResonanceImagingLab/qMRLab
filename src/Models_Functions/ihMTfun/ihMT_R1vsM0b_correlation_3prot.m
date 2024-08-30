function ihMT_R1vsM0b_correlation_3prot(obj)
%% Generate the M0B mapping to R1 from simulation results and acquired data

%  NOTE: niak_read_vol replaced with minc_read  %
%  NOTE: niak_write_vol replaced with minc_write %

%% Load images:

% Where Images are you will be using 
DATADIR = obj.options.R1vsM0bMapping_DataDirectory;
%DATADIR = '/Users/amiedemmans/Desktop/GitHub/Image_files/ihMT_images_from_chris/matlab/';


% Where outputs will be saved/ previous results are 
OutputDir = obj.options.R1vsM0bMapping_OutputDirectory;
%OutputDir = '/Users/amiedemmans/Desktop/GitHub/qMRLab_test_dir/Test1';


%image names:
%mtw_fn = {'dual_reg.mnc' 'pos_reg.mnc' 'neg_reg.mnc' 'sparseMP2RAGE_M0.mnc.gz' 'sparseMP2RAGE_T1.mnc.gz'};
mtw_fn = {data.MTw_dual, data.MTw_single_pos, data.MTw_single_neg, data.T1map, data.S0map};

for i = 1:size(mtw_fn,2)
    fn = fullfile(DATADIR,mtw_fn{i});
    [hdr, img] = minc_read(fn); % niak_read_vol
    comb_mtw(:,:,:,i) = img; %.img;
end


%% Load the mask
fn = fullfile(DATADIR,'itkmask.mnc');
[~, mask] = minc_read(fn);
%mask1 = permute(mask,[2 3 1]); % conversion between minc and nii reorients it

% If making a mask 
% mask = zeros(size(lfa)); 
% mask(lfa>175) = 1;

%% Some B1 issues so lets try and load that
[~, b1] = minc_read(fullfile(DATADIR,'resampled_b1field.mnc')); 
b1 = double(b1);
b1 = limitHandler(b1, 0.5,1.4);
b1 = permute(b1, [3 1 2]);
b1 = ihMT_imgaussfilt3_withMask(b1, mask, 5); %light smoothing to the map

%figure; imshow3Dfull(b1, [0.6 1.2],jet)

%% List the images 
dual = comb_mtw(:,:,:,1);
pos = comb_mtw(:,:,:,2);
neg = comb_mtw(:,:,:,3);

spApp_mp2 = comb_mtw(:,:,:,4); 
spT1_map = comb_mtw(:,:,:,5); 


%% Mask -> bet result touched up in itk, then threshold CSF and some dura
% AMIE- move to just before M0b maps  
%mask = mask1;
% mask(spT1_map > 2500) = 0;
% mask(spT1_map < 650) = 0;
% mask(isnan(spT1_map)) = 0;
% mask = bwareaopen(mask, 10000,6);
% figure; imshow3Dfullseg(spT1_map, [300 2500],mask)


%% Compute MTsat ihMTsat

%% Protocol 
% flipA = obj.Prot.PulseSequenceParams.Mat(3);
% TR = obj.Prot.PulseSequenceParams.Mat(4)*1000;
% DummyEcho = obj.Prot.PulseSequenceParams.Mat(10);
% echoSpacing = obj.Prot.PulseSequenceParams.Mat(14)*1000;
% numExcitation = obj.Prot.PulseSequenceParams.Mat(6) + DummyEcho;
flipA = 6; 
TR = 1.14*1000;
%TR = 23;
DummyEcho = 2; 
echoSpacing = 0.0077*1000;
numExcitation = 10-DummyEcho;

sat_dual = ihMT_calcMTsatThruLookupTablewithDummyV3( dual, b1, spT1_map, mask,spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_pos  = ihMT_calcMTsatThruLookupTablewithDummyV3( pos, b1, spT1_map, mask, spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_neg  = ihMT_calcMTsatThruLookupTablewithDummyV3( neg, b1, spT1_map, mask, spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);

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
%% This isn't working because dividing by zero --> messing everything up %%
% getting values of inf so when get to line 164 it doesnt work
%R1_s = (1./spT1_map);
R1_s(isinf(R1_s)) = 0; % This I think fixed it 

% initialize matrices
M0b_dual = zeros(size(sat_dual));
M0b_pos = zeros(size(sat_dual));
M0b_neg = zeros(size(sat_dual));


%% SPEED IT UP BY DOING A FEW AXIAL SLICES
axialStart = 126; % 65
axialStop = axialStart+3;%115;
% check
figure; imshow3Dfull(sat_dual(:,axialStart:axialStop,:) , [0 0.06], jet)

%b1_1 = obj.Prot.PulseSequenceParams.Mat(8);
b1_1 = 11.6; %AMIE
tic %  ~ 2hrs to run 
for i = 1:size(sat_dual,1) % went to 149
    
    for j = axialStart:axialStop % 1:size(sat_dual,2) % for axial slices
        for k =  1:size(sat_dual,3) % sagital slices  65
            
            if mask(i,j,k) > 0 %&& dual_s(i,j,k,3) > 0
                                
                 [M0b_dual(i,j,k), ~,  ~]  = ihMT_fit_M0b_v2( b1(i,j,k), R1_s(i,j,k), sat_dual(i,j,k), fitValues_D);
                 [M0b_pos(i,j,k),  ~,  ~]  = ihMT_fit_M0b_v2( b1_1*b1(i,j,k), R1_s(i,j,k), sat_pos(i,j,k), fitValues_S);               
                 [M0b_neg(i,j,k),  ~,  ~]  = ihMT_fit_M0b_v2( b1_1*b1(i,j,k), R1_s(i,j,k), sat_neg(i,j,k), fitValues_S);
                 
            end
        end
    end
    disp(i/size(sat_dual,1) *100)
    toc 
end
[~, M0b_dual] = minc_read('M0b_dual.mnc.gz');
[~, M0b_pos] = minc_read('M0b_pos.mnc.gz');
[~, M0b_neg] = minc_read('M0b_neg.mnc.gz');

%% this took 30hours for 1mm isotropic full brain dataset. * was running fitting in another matlab
    % instance, so could be easily sped up running on its own and/or adding
    % the parfor loop. 


figure; imshow3Dfull(M0b_pos, [0 0.15],jet)
figure; imshow3Dfull(M0b_neg, [0 0.15], jet)  
figure; imshow3Dfull(M0b_dual, [0 0.15],jet)

% export
mkdir(fullfile(OutputDir,'processing'))
hdr.file_name = fullfile(OutputDir,'processing/M0b_dual.mnc.gz'); 
minc_write(hdr.file_name, hdr, M0b_dual); % Need to add filename for each minc_write %
hdr.file_name = fullfile(OutputDir,'processing/M0b_pos.mnc.gz');
minc_write(hdr.file_name, hdr, M0b_pos);
hdr.file_name = fullfile(OutputDir,'processing/M0b_neg.mnc.gz'); 
minc_write(hdr.file_name, hdr, M0b_neg);


%% With M0B maps made, correlate with R1 and update the fitValues file. 

% use this fake mask to get rid of dura. 
tempMask = mask;
tempMask = imerode(tempMask, strel('sphere',2));
figure; imshow3Dfullseg(M0b_dual, [0 0.15],tempMask)

mkdir(fullfile(OutputDir,'figures'));

% Optimized Approach
fitValues_D  = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_dual, tempMask, fitValues_D, fullfile(OutputDir,'figures/R1vsM0b_dual.png'), fullfile(OutputDir,'fitValues_D.mat'));
fitValues_SP = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_pos, tempMask, fitValues_S, fullfile(OutputDir,'figures/R1vsM0b_pos.png'), fullfile(OutputDir,'fitValues_SP.mat'));
fitValues_SN = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_neg, tempMask, fitValues_S, fullfile(OutputDir,'figures/R1vsM0b_neg.png'), fullfile(OutputDir,'fitValues_SN.mat'));


%% Now use these results to B1 correct the data:
OutputDir = DATADIR;

b1_1 = obj.Prot.PulseSequenceParams.Mat(8);

corr_prot_d = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_D);
corr_prot_p = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_SP);
corr_prot_n = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_SN);


% Part 2, apply correction map
sat_dual_c = (sat_dual + sat_dual.* corr_prot_d) .* mask1;
sat_pos_c  = (sat_pos + sat_pos.* corr_prot_p) .* mask1;
sat_neg_c  = (sat_neg + sat_neg.* corr_prot_n) .* mask1;
ihmt_c      = sat_dual_c - (sat_pos_c + sat_neg_c)/2;


ihmt_c = double(limitHandler(ihmt_c,0, 0.05));

ihmt_c( ihmt_c >= 0.05) = 0;

%% View results
figure; imshow3Dfull(sat_dual1_c , [0 0.06], jet); 
figure; imshow3Dfull(sat_pos1_c , [0 0.06], jet)
figure; imshow3Dfull(sat_neg1_c , [0 0.06], jet); 

figure; imshow3Dfull(ihmt_c , [0 0.06], jet)

%% Other things, save if you want
mkdir(strcat(DATADIR,'R1vsM0b_results') )
 
hdr.file_name = strcat(DATADIR,'R1vsM0b_results/MTsat_dual_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_dual_c);
hdr.file_name = strcat(DATADIR,'R1vsM0b_results/MTsat_pos_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_pos_c);
hdr.file_name = strcat(DATADIR,'R1vsM0b_results/MTsat_neg_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_neg_c);
hdr.file_name = strcat(DATADIR,'R1vsM0b_results/ihMTsat_1.mnc.gz'); minc_write(hdr.file_name, hdr, ihmt_c);

hdr.file_name = strcat(DATADIR,'R1vsM0b_results/b1.mnc.gz'); minc_write(hdr.file_name, hdr, b1);

%% 

ihmtSlice1 = squeeze( ihmt_c(:,126,:));

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










