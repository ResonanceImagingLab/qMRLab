%% Correct MTsat maps from 3 protocols after running:
%   simSeq_M0B_R1obs_3prot.m  and...
%   CR_R1vsM0B_correlation.m
function mtsat_c = mt_sat_b1sim_correctMTsat(data, fitValues_SP, flipA, TR, DummyEcho, echoSpacing, numExcitation)

%% Load images:
pos = data.pos; 

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

mtsat  = ihMT_calcMTsatThruLookupTablewithDummyV3( pos, b1, T1_map, mask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% figure; imshow3Dfull(sat_dual , [0 0.06], jet); figure; imshow3Dfull(sat_pos1 , [0 0.05], jet);  figure; imshow3Dfull(sat_neg1 , [0 0.05], jet); 

%% Now use these results to B1 correct the data:
%OutputDir = DATADIR;
if isfield(fitValues_SP, 'fitValues')
    fitValues_SP = fitValues_SP.fitValues;
end

R1_s = (1./T1_map) *1000; 
R1_s(isinf(R1_s)) = 0; 

corr_p = MTsat_B1corr_factor_map(b1, R1_s, 1, fitValues_SP);

% Part 2, apply correction map

mtsat_c  = (mtsat + mtsat.* corr_p) .* mask;


mtsat_c = double(limitHandler(mtsat_c, 0, 0.3));

%% View results

% figure; imshow3Dfull(mtsat_c , [0 0.06], jet)























