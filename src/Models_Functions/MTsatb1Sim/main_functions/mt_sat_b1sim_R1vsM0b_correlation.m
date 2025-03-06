function fitValues_SP = mt_sat_b1sim_R1vsM0b_correlation(data, fitValues_single, flipA, TR, DummyEcho, echoSpacing, numExcitation, OutputDir, scaleS0)

%% Load the images 
mtw = data.MTw; 
S0_map = data.S0map;
T1_map = data.T1map;
b1 = data.b1; 

if scaleS0
    S0_map = S0_map/2.5;
    disp('dividing S0 map by 2.5 to account for scanner gain differences between GRE and MP2RAGE (Siemens ONLY)')
end

% Load the mask
if ~isempty(data.mask)
    mask = data.mask; 
    maskFlag = 0;
else 
    mask = zeros(size(mtw));
    mask(b1>0) = 1; 
    maskFlag = 1;
end 

%% Compute ihMTsat 

tempMask = mask;
tempMask(T1_map > 2500) = 0;
tempMask(T1_map < 650) = 0;
tempMask(isnan(T1_map)) = 0;

sat_mtw  = ihMT_calcMTsatThruLookupTablewithDummyV3( mtw, b1, T1_map, tempMask, S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);

figure; imshow3Dfull(sat_mtw , [0 0.05], jet);
disp('Check to make sure the MTsat map values look appropriate. If not, check scaling of S_0. If they look good, press any key to continue...')
pause;


% load in the fit results for VFA - Optimal
fitValues_single = fitValues_single.fitValues;

R1_s = (1./T1_map) *1000; 
R1_s(isinf(R1_s)) = 0;  

% initialize matrices
M0b_mtw = zeros(size(sat_mtw));

disp('To speed things up, we will only run the fitting on a middle block of slices');
zSlice = round(size(sat_mtw,3)/2);

ww = 0; % Waitbar counter
h = [];
stop = 0;
% Create waitbar
if (~exist('wait','var') || isempty(wait))
    wait = 1;   % waitbar is on by default
end

if (wait)
    h = waitbar(0,'','Name','Simulating data','CreateCancelBtn',...
        'if ~strcmp(get(gcbf,''Name''),''canceling...''), setappdata(gcbf,''canceling'',1); set(gcbf,''Name'',''canceling...''); else delete(gcbf); end');
    setappdata(h,'canceling',0)
    setappdata(0,'Cancel',0);
end

denom = size(sat_mtw,1);
tic %  ~ 2hrs to run 
for i = 1:denom 

    if (wait)
        % Update waitbar
        ww = ww+1;
        waitbar(ww/size(sat_mtw,1), h, sprintf(' %.1f percent complete', i/denom*100) );
    end
    
    for j = 1:size(sat_mtw,2)  % for axial slices
        for k =  zSlice-10:zSlice+10 % sagital slices  
            
            if tempMask(i,j,k) > 0
                                
                 [M0b_mtw(i,j,k),  ~]  = ihMT_fit_M0b_v2( b1(i,j,k), R1_s(i,j,k), sat_mtw(i,j,k), fitValues_single);               
                 
            end
            if (wait)
                % Allows user to cancel
                if getappdata(h,'canceling')
                    stop = 1;
                    setappdata(0,'Cancel',1);
                    break;
                end
            end
        end

        if (stop)
                delete(h);
                error('Simulations Cancelled');
        end

    end 
end
disp('Fitting Finished. Total time: ')
toc
disp('Now Calculating Values...')
% Delete waitbar
delete(h);


%figure('WindowStyle', 'docked'); imshow3Dfull(M0b_mtw, [0 0.15],jet)

% export
% mkdir(fullfile(OutputDir,'processing'))
% hdr.file_name = fullfile(OutputDir,'processing/M0b_mtw.mnc.gz');
% minc_write(hdr.file_name, hdr, M0b_mtw);

%% With M0B maps made, correlate with R1 and update the fitValues file. 

% If not using own mask, increase imerode to sphere 4 
tempMask = bwareaopen(tempMask, 10000,6);
if maskFlag 
    tempMask = imerode(tempMask, strel('sphere',4));
else 
    tempMask = imerode(tempMask, strel('sphere',2));
end

mkdir(fullfile(OutputDir,'figures'));
mkdir(fullfile(OutputDir, 'R1vsM0b_results'))

% Optimized Approach
fitValues_SP = ihMT_generate_R1vsM0B_correlation( R1_s, M0b_mtw, tempMask, fitValues_single, fullfile(OutputDir,'figures/R1vsM0b_MTsat.png'), fullfile(OutputDir,'R1vsM0b_results/fitValues_SP.mat'));


















