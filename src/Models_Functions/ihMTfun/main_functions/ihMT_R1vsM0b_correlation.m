function [fitValues_Dual, fitValues_SP, fitValues_SN] = ihMT_R1vsM0b_correlation(data, fitValues_dual, fitValues_single, flipA, TR, DummyEcho, echoSpacing, numExcitation, OutputDir)

%% Load the images 
dual = data.dual; 
pos = data.pos; 
neg = data.neg; 

S0_map = data.M0map;
T1_map = data.T1map;
b1 = data.b1; 

% Load the map 
if ~isempty(data.mask)
    mask = data.mask; 
    maskFlag = 0;
else 
    mask = zeros(size(dual));
    mask(b1>0) = 1; 
    maskFlag = 1;
end 

%% Compute ihMTsat 

tempMask = mask;
tempMask(T1_map > 2500) = 0;
tempMask(T1_map < 650) = 0;
tempMask(isnan(T1_map)) = 0;

sat_dual = ihMT_calcMTsatThruLookupTablewithDummyV3( dual, b1, T1_map, tempMask,S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_pos  = ihMT_calcMTsatThruLookupTablewithDummyV3( pos, b1, T1_map, tempMask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_neg  = ihMT_calcMTsatThruLookupTablewithDummyV3( neg, b1, T1_map, tempMask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% figure; imshow3Dfull(sat_dual , [0 0.06], jet); figure; imshow3Dfull(sat_pos1 , [0 0.05], jet);  figure; imshow3Dfull(sat_neg1 , [0 0.05], jet); 

% load in the fit results for VFA - Optimal
fitValues_single = fitValues_single.fitValues;
fitValues_dual = fitValues_dual.fitValues;

R1_s = (1./T1_map) *1000; 
R1_s(isinf(R1_s)) = 0;  

% initialize matrices
M0b_dual = zeros(size(sat_dual));
M0b_pos = zeros(size(sat_dual));
M0b_neg = zeros(size(sat_dual));

disp('Code will take ~ 3 hours to run');
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

%figure('WindowStyle', 'docked'); imshow3Dfull(M0b_pos, [0 0.15],jet)
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

% If not using own mask, increase imerode to sphere 4 
tempMask = bwareaopen(tempMask, 10000,6);
if maskFlag 
    tempMask = imerode(tempMask, strel('sphere',4));
else 
    tempMask = imerode(tempMask, strel('sphere',2));
end


% figure; imshow3Dfullseg(M0b_dual, [0 0.15],tempMask)

mkdir(fullfile(OutputDir,'figures'));
mkdir(fullfile(OutputDir, 'R1vsM0b_results'))


% Optimized Approach
fitValues_Dual  = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_dual, tempMask, fitValues_dual, fullfile(OutputDir,'figures/R1vsM0b_dual.png'), fullfile(OutputDir,'R1vsM0b_results/fitValues_Dual.mat'));
fitValues_SP = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_pos, tempMask, fitValues_single, fullfile(OutputDir,'figures/R1vsM0b_pos.png'), fullfile(OutputDir,'R1vsM0b_results/fitValues_SP.mat'));
fitValues_SN = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_neg, tempMask, fitValues_single, fullfile(OutputDir,'figures/R1vsM0b_neg.png'), fullfile(OutputDir,'R1vsM0b_results/fitValues_SN.mat'));


















