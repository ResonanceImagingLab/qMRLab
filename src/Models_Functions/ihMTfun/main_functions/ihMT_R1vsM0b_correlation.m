function [fitValues_Dual, fitValues_SP, fitValues_SN] = ihMT_R1vsM0b_correlation(obj, data, fitValues_dual, fitValues_single)

% Sequence Simulations Results Directory 
%SeqSimDir = obj.options.R1vsM0bMapping_DataDirectory;

% Directory where results will be saved 
OutputDir =  obj.options.R1vsM0bMapping_SeqSimDirectory;

%% Load the images 
dual = data.dual; 
pos = data.pos; 
neg = data.neg; 

S0_map = data.M0map;
T1_map = data.T1map;
b1 = data.b1; 

% Load the map 
if isfield(data, 'mask')
    mask = data.mask; 
else 
    mask = zeros(size(dual));
    mask(b1>0) = 1; 
end 


%% Protocol
flipA = obj.Prot.PulseSequenceParams.Mat(3);
TR = obj.Prot.PulseSequenceParams.Mat(4); % ms
DummyEcho = obj.Prot.PulseSequenceParams.Mat(10);
echoSpacing = obj.Prot.PulseSequenceParams.Mat(14); % ms 
numExcitation = obj.Prot.PulseSequenceParams.Mat(6) + DummyEcho;

%% Compute ihMTsat 

tempMask = mask;
tempMask(T1_map > 2500) = 0;
tempMask(T1_map < 650) = 0;
tempMask(isnan(T1_map)) = 0;

sat_dual = ihMT_calcMTsatThruLookupTablewithDummyV3( dual, b1, T1_map, tempMask,S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_pos  = ihMT_calcMTsatThruLookupTablewithDummyV3( pos, b1, T1_map, tempMask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_neg  = ihMT_calcMTsatThruLookupTablewithDummyV3( neg, b1, T1_map, tempMask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% figure; imshow3Dfull(sat_dual , [0 0.06], jet); figure; imshow3Dfull(sat_pos , [0 0.05], jet);  figure; imshow3Dfull(sat_neg , [0 0.05], jet); 

% load in the fit results for VFA - Optimal
fitValues_single = fitValues_single.fitValues;
fitValues_dual = fitValues_dual.fitValues;

R1_s = (1./T1_map) *1000; 
R1_s(isinf(R1_s)) = 0;  

% initialize matrices
M0b_dual = zeros(size(sat_dual));
M0b_pos = zeros(size(sat_dual));
M0b_neg = zeros(size(sat_dual));

% Speed up by doing only a few axial slices 
% axialStart = 126; % 65
% axialStop = axialStart+3;%115;
%figure; imshow3Dfull(sat_dual(:,axialStart:axialStop,:) , [0 0.06], jet)

%disp('Code will take ~ 2 hours to run');
tic %  ~ 2hrs to run 
for i = 1:size(sat_dual,1) % went to 149
    
    for j = 1:size(sat_dual,2) %j = axialStart:axialStop  % % for axial slices
        for k =  1:size(sat_dual,3) % sagital slices  65
            
            if tempMask(i,j,k) > 0 %&& dual_s(i,j,k,3) > 0
                                
                 [M0b_dual(i,j,k), ~]  = ihMT_fit_M0b_v2( b1(i,j,k), R1_s(i,j,k), sat_dual(i,j,k), fitValues_dual);
                 [M0b_pos(i,j,k),  ~]  = ihMT_fit_M0b_v2( b1(i,j,k), R1_s(i,j,k), sat_pos(i,j,k), fitValues_single);               
                 [M0b_neg(i,j,k),  ~]  = ihMT_fit_M0b_v2( b1(i,j,k), R1_s(i,j,k), sat_neg(i,j,k), fitValues_single);
                 
            end
        end
    end
    disp(i/size(sat_dual,1) *100)
    toc 
end

%figure('WindowStyle', 'docked') % docked in matlab i believe 

% figure('WindowStyle', 'docked'); imshow3Dfull(M0b_pos, [0 0.15],jet)
% figure('WindowStyle', 'docked'); imshow3Dfull(M0b_neg, [0 0.15], jet)  
% figure('WindowStyle', 'docked'); imshow3Dfull(M0b_dual, [0 0.15],jet)

% export
% mkdir(fullfile(OutputDir,'processing'))
% hdr.file_name = fullfile(OutputDir,'processing/M0b_dual.mnc.gz'); 
% minc_write(hdr.file_name, hdr, M0b_dual); % Need to add filename for each minc_write %
% hdr.file_name = fullfile(OutputDir,'processing/M0b_pos.mnc.gz');
% minc_write(hdr.file_name, hdr, M0b_pos);
% hdr.file_name = fullfile(OutputDir,'processing/M0b_neg.mnc.gz'); 
% minc_write(hdr.file_name, hdr, M0b_neg);

%% With M0B maps made, correlate with R1 and update the fitValues file. 

% use this fake mask to get rid of dura. 
% tempMask = mask;
% tempMask(T1_map > 2500) = 0;
% tempMask(T1_map < 650) = 0;
% tempMask(isnan(T1_map)) = 0;
tempMask = bwareaopen(tempMask, 10000,6);
tempMask = imerode(tempMask, strel('sphere',2));
% figure; imshow3Dfullseg(M0b_dual, [0 0.15],tempMask)

mkdir(fullfile(OutputDir,'figures2'));
mkdir(fullfile(OutputDir, 'R1vsM0b_results2'))

%OutputDir = '/Users/amiedemmans/Documents/ihMT_Tests/Test4/';

% Optimized Approach
fitValues_Dual  = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_dual, tempMask, fitValues_dual, fullfile(OutputDir,'figures2/R1vsM0b_dual.png'), fullfile(OutputDir,'R1vsM0b_results2/fitValues_Dual.mat'));
fitValues_SP = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_pos, tempMask, fitValues_single, fullfile(OutputDir,'figures/R1vsM0b_pos.png'), fullfile(OutputDir,'R1vsM0b_results/fitValues_SP.mat'));
fitValues_SN = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_neg, tempMask, fitValues_single, fullfile(OutputDir,'figures/R1vsM0b_neg.png'), fullfile(OutputDir,'R1vsM0b_results/fitValues_SN.mat'));

%% This is to be just placed in MTsat code 
% --> will need to add if statements probably for R1_s
%% Now use these results to B1 correct the data:
%OutputDir = DATADIR;

% b1_1 = obj.Prot.PulseSequenceParams.Mat(8);
% 
% corr_prot_d = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_D);
% corr_prot_p = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_SP);
% corr_prot_n = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_SN);
% 
% 
% % Part 2, apply correction map
% sat_dual_c = (sat_dual + sat_dual.* corr_prot_d) .* mask1;
% sat_pos_c  = (sat_pos + sat_pos.* corr_prot_p) .* mask1;
% sat_neg_c  = (sat_neg + sat_neg.* corr_prot_n) .* mask1;
% ihmt_c      = sat_dual_c - (sat_pos_c + sat_neg_c)/2;
% 
% 
% ihmt_c = double(limitHandler(ihmt_c,0, 0.05));
% 
% ihmt_c( ihmt_c >= 0.05) = 0;
% 
% %% View results
% figure; imshow3Dfull(sat_dual1_c , [0 0.06], jet); 
% figure; imshow3Dfull(sat_pos1_c , [0 0.06], jet)
% figure; imshow3Dfull(sat_neg1_c , [0 0.06], jet); 
% 
% figure; imshow3Dfull(ihmt_c , [0 0.06], jet)
% 
% %% Other things, save if you want
% mkdir(strcat(OutputDir,'R1vsM0b_results') )
% 
% hdr.file_name = strcat(OutputDir,'R1vsM0b_results/MTsat_dual_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_dual_c);
% hdr.file_name = strcat(OutputDir,'R1vsM0b_results/MTsat_pos_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_pos_c);
% hdr.file_name = strcat(OutputDir,'R1vsM0b_results/MTsat_neg_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_neg_c);
% hdr.file_name = strcat(OutputDir,'R1vsM0b_results/ihMTsat_1.mnc.gz'); minc_write(hdr.file_name, hdr, ihmt_c);
% 
% hdr.file_name = strcat(OutputDir,'R1vsM0b_results/b1.mnc.gz'); minc_write(hdr.file_name, hdr, b1);
% 
% %% 
% 
% ihmtSlice1 = squeeze( ihmt_c(:,126,:));
% 
% figure; imagesc(ihmtSlice1); axis image;
% colormap(gray)
% caxis([0 0.03])
% hold on
% line([0,175], [165,166], 'Color', 'r');
% 
% 
% ihmtProf1 = ihmtSlice1(165,:);
% 
% % normalize:
% ihmtProf1 = ihmtProf1/ max(ihmtProf1);
% 
% figure
% plot(ihmtProf1,'LineWidth',2); 
% title('Line Profile (L-R)');
% xlim([15, 160])
% hold off
% xlabel('Voxel Index');
% ylabel('Relative ihMT_{sat}')
% ax = gca; ax.FontSize = 20; 
% set(gcf,'position',[10,400,1000,600])   
% 

















