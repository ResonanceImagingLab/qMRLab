%% Correct MTsat maps from 3 protocols after running:
%   simSeq_M0B_R1obs_3prot.m  and...
%   CR_R1vsM0B_correlation.m
function [sat_dual_c, sat_pos_c, sat_neg_c, ihmt_c]=ihMT_correctMTsat_3prot(data, fitValues_D, fitValues_SP, fitValues_SN, flipA, TR, DummyEcho, echoSpacing, numExcitation)

%% Load images:

dual = data.dual; 
pos = data.pos; 
neg = data.neg; 

S0_map = data.M0map;
T1_map = data.T1map;
b1 = data.b1; 

% Load the map 
if ~isempty(data.mask)
    mask = data.mask; 

else 
    mask = zeros(size(dual));
    mask(b1>0) = 1; 

end 

%% Mask -> bet result touched up in itk, then threshold CSF and some dura

mask(T1_map > 2500) = 0;
mask(T1_map < 500) = 0;
mask(isnan(T1_map)) = 0;

%% Compute MTsat ihMTsat

sat_dual = ihMT_calcMTsatThruLookupTablewithDummyV3( dual, b1, T1_map, mask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_pos  = ihMT_calcMTsatThruLookupTablewithDummyV3( pos, b1, T1_map, mask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_neg  = ihMT_calcMTsatThruLookupTablewithDummyV3( neg, b1, T1_map, mask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% figure; imshow3Dfull(sat_dual , [0 0.06], jet); figure; imshow3Dfull(sat_pos1 , [0 0.05], jet);  figure; imshow3Dfull(sat_neg1 , [0 0.05], jet); 

%% Now use these results to B1 correct the data:
%OutputDir = DATADIR;
fitValues_D = fitValues_D.fitValues;
fitValues_SP = fitValues_SP.fitValues;
fitValues_SN = fitValues_SN.fitValues;

R1_s = (1./T1_map) *1000; 
R1_s(isinf(R1_s)) = 0; 

corr_d = MTsat_B1corr_factor_map(b1, R1_s, 1, fitValues_D);
corr_p = MTsat_B1corr_factor_map(b1, R1_s, 1, fitValues_SP);
corr_n = MTsat_B1corr_factor_map(b1, R1_s, 1, fitValues_SN);

% Part 2, apply correction map
sat_dual_c = (sat_dual + sat_dual.* corr_d) .* mask;
sat_pos_c  = (sat_pos + sat_pos.* corr_p) .* mask;
sat_neg_c  = (sat_neg + sat_neg.* corr_n) .* mask;
ihmt_c      = sat_dual_c - (sat_pos_c + sat_neg_c)/2;

ihmt_c = double(limitHandler(ihmt_c,0, 0.15));

%% View results

% figure; imshow3Dfull(sat_dual_c , [0 0.06], jet); 
% figure; imshow3Dfull(sat_pos_c , [0 0.06], jet)
% figure; imshow3Dfull(sat_neg_c , [0 0.06], jet); 
% 
% figure; imshow3Dfull(ihmt_c , [0 0.02], jet);

%% Other things, save if you want

% hdr.file_name = strcat(DATADIR,'matlab/MTsat_dual_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_dual_c);
% hdr.file_name = strcat(DATADIR,'matlab/MTsat_pos_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_pos_c);
% hdr.file_name = strcat(DATADIR,'matlab/MTsat_neg_1.mnc.gz'); minc_write(hdr.file_name, hdr, sat_neg_c);
% hdr.file_name = strcat(DATADIR,'matlab/ihMTsat_1.mnc.gz'); minc_write(hdr.file_name, hdr, ihmt_c);
% 
% hdr.file_name = strcat(DATADIR,'matlab/b1.mnc.gz'); minc_write(hdr.file_name, hdr, b1);

ihmtSlice1 = squeeze( ihmt_c(:,126,:));

figure; imagesc(ihmtSlice1); axis image;
colormap(gray)
clim([0 0.03]) % was caxis
hold on
line([0,175], [165,166], 'Color', 'r');


ihmtProf1 = ihmtSlice1(165,:); % These are all zeros 

% normalize:
ihmtProf1 = ihmtProf1/ max(ihmtProf1);  % This produces Nan 

figure
plot(ihmtProf1,'LineWidth',2); 
title('Line Profile (L-R)');
%xlim([15, 160])
hold off
xlabel('Voxel Index');
ylabel('Relative ihMT_{sat}')
ax = gca; ax.FontSize = 20; 
set(gcf,'position',[10,400,1000,600])

























