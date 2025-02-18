function fitValues_SP = mt_sat_b1sim_R1vsM0b_correlation(data, fitValues_dual, fitValues_single, flipA, TR, DummyEcho, echoSpacing, numExcitation, OutputDir)

%% Load the images 
pos = data.pos; 
S0_map = data.M0map;
T1_map = data.T1map;
b1 = data.b1; 

% Load the mask
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

sat_pos  = ihMT_calcMTsatThruLookupTablewithDummyV3( pos, b1, T1_map, tempMask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% figure; imshow3Dfull(sat_pos , [0 0.04], jet);

% load in the fit results for VFA - Optimal
fitValues_single = fitValues_single.fitValues;

R1_s = (1./T1_map) *1000; 
R1_s(isinf(R1_s)) = 0;  

% initialize matrices
M0b_pos = zeros(size(sat_pos));

disp('Code will take ~ 3 hours to run');
tic %  ~ 2hrs to run 
for i = 1:size(sat_pos,1) 
    
    for j = 1:size(sat_pos,2)  % for axial slices
        for k =  1:size(sat_pos,3) % sagital slices  
            
            if tempMask(i,j,k) > 0 %&& dual_s(i,j,k,3) > 0
                                
                 [M0b_pos(i,j,k),  ~]  = ihMT_fit_M0b_v2( b1(i,j,k), R1_s(i,j,k), sat_pos(i,j,k), fitValues_single);               
                 
            end
        end
    end
    disp(i/size(sat_pos,1) *100)
    toc 
end

%figure('WindowStyle', 'docked'); imshow3Dfull(M0b_pos, [0 0.15],jet)

% export
% mkdir(fullfile(OutputDir,'processing'))
% hdr.file_name = fullfile(OutputDir,'processing/M0b_pos.mnc.gz');
% minc_write(hdr.file_name, hdr, M0b_pos);

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
fitValues_SP = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_pos, tempMask, fitValues_single, fullfile(OutputDir,'figures/R1vsM0b_pos.png'), fullfile(OutputDir,'R1vsM0b_results/fitValues_SP.mat'));


















