%% Correct MTsat maps from 3 protocols after running:
%   simSeq_M0B_R1obs_3prot.m  and...
%   CR_R1vsM0B_correlation.m
function [sat_dual_c, sat_pos_c, sat_neg_c, ihmt_c]=ihMT_correctMTsat_3prot(obj,data)

% Data Directory 
DATADIR = obj.options.RunihMTsatCalculation_DataDirectory;

% Directory where results will be saved 
OutputDir = obj.options.RunihMTsatCalculation_OutputDirectory;

%% Load images:

dual = data.dual_reg; 
pos = data.pos_reg; 
neg = data.neg_reg; 

S0_map = data.sparseMP2RAGE_M0;
T1_map = data.sparseMP2RAGE_T1;
b1 = data.b1_permute; 

% Load the map 
if isfield(data, 'mask')
    mask = data.mask; 
else 
    mask = zeros(size(dual));
    mask(b1>0) = 1; 
end 

%% Mask -> bet result touched up in itk, then threshold CSF and some dura

% mask = mask1;
% mask(spT1_map > 2500) = 0;
% mask(spT1_map < 500) = 0;
% mask(isnan(spT1_map)) = 0;
% mask = bwareaopen(mask, 10000,6);
% figure; imshow3Dfullseg(spT1_map, [300 2500],mask1)

%% Compute MTsat ihMTsat

%% Protocol 1
flipA = obj.Prot.PulseSequenceParams.Mat(3);
TR = obj.Prot.PulseSequenceParams.Mat(4)*1000;
DummyEcho = obj.Prot.PulseSequenceParams.Mat(10);
echoSpacing = obj.Prot.PulseSequenceParams.Mat(14)*1000;
numExcitation = obj.Prot.PulseSequenceParams.Mat(6) + DummyEcho;

sat_dual = ihMT_calcMTsatThruLookupTablewithDummyV3( dual, b1, T1_map, mask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_pos  = ihMT_calcMTsatThruLookupTablewithDummyV3( pos, b1, T1_map, mask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_neg  = ihMT_calcMTsatThruLookupTablewithDummyV3( neg, b1, T1_map, mask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% figure; imshow3Dfull(sat_dual1 , [0 0.06], jet); figure; imshow3Dfull(sat_pos1 , [0 0.05], jet);  figure; imshow3Dfull(sat_neg1 , [0 0.05], jet); 

%% load in the fit results
fitValues_D = load(fullfile(OutputDir,'fitValues_D.mat'));
fitValues_D = fitValues_D.fitValues;
fitValues_SP = load(fullfile(OutputDir,'fitValues_SP.mat'));
fitValues_SP = fitValues_SP.fitValues;
fitValues_SN = load(fullfile(OutputDir,'fitValues_SN.mat'));
fitValues_SN = fitValues_SN.fitValues;


%% Now use these results to B1 correct the data:
%OutputDir = DATADIR;

b1_1 = 11.4;

corr_d = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_D);
corr_p = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_SP);
corr_n = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_SN);

% Part 2, apply correction map
sat_dual_c = (sat_dual + sat_dual.* corr_d) .* mask;
sat_pos_c  = (sat_pos + sat_pos.* corr_p) .* mask;
sat_neg_c  = (sat_neg + sat_neg.* corr_n) .* mask;
ihmt_c      = sat_dual_c - (sat_pos_c + sat_neg_c)/2;

ihmt_c = double(limitHandler(ihmt_c,0, 0.15));

%% View results
% Can maybe get rid of this assuming corrected image will appear in qMRLab 

figure; imshow3Dfull(sat_dual_c , [0 0.06], jet); 
figure; imshow3Dfull(sat_pos_c , [0 0.06], jet)
figure; imshow3Dfull(sat_neg_c , [0 0.06], jet); 

figure; imshow3Dfull(ihmt_c , [0 0.02], jet);

%% Other things, save if you want

hdr.file_name = strcat(DATADIR,'matlab/MTsat_dual_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_dual_c);
hdr.file_name = strcat(DATADIR,'matlab/MTsat_pos_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_pos_c);
hdr.file_name = strcat(DATADIR,'matlab/MTsat_neg_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_neg_c);
hdr.file_name = strcat(DATADIR,'matlab/ihMTsat_1.mnc.gz'); minc_write(hdr.file_name, hdr, ihmt_c);

hdr.file_name = strcat(DATADIR,'matlab/b1.mnc.gz'); minc_write(hdr.file_name, hdr, b1);

ihmtSlice1 = squeeze( ihmt_c(:,126,:));

figure; imagesc(ihmtSlice1); axis image;
colormap(gray)
clim([0 0.03]) % was caxis
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
























