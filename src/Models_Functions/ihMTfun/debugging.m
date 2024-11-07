%% Test script for ihMT_fit_M0b_v2.m 

B1_ref = 0.9520; % disp value 

msat = 0.0342; % disp value 

flipA = 6; 
TR = 100;
DummyEcho = 2; 
echoSpacing = 7.66; 
numExcitation = 10;

fitValues_dual = load('/Users/amiedemmans/Documents/ihMT_Tests/Test2/fitValues_D.mat');
fitValues = fitValues_dual.fitValues;

[~, dual] = minc_read('dual_reg.mnc'); 
[~, T1_map] = minc_read('T1map.mnc');
[~, S0_map] = minc_read('M0map.mnc'); 
[~, mask] = minc_read('mask.mnc'); 
[~, B1_ref] = minc_read('b1_file.mnc');


tempMask = mask;
tempMask(T1_map > 2500) = 0;
tempMask(T1_map < 650) = 0;
tempMask(isnan(T1_map)) = 0;

Raobs = (1./T1_map) *1000; 
Raobs(isinf(Raobs)) = 0;  


sat_dual = ihMT_calcMTsatThruLookupTablewithDummyV3( dual, B1_ref, T1_map, tempMask,S0_map, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% i = 1:size(sat_dual, 1);
% j = 126:129; 
% k = 1:size(sat_dual, 3);


% i = 126:129;
% j = 126:129; 
% k = 126:129;


i = 64; j = 66; k = 46;

B1_ref = B1_ref(i, j, k, :); 
msat = sat_dual(i, j, k, :);
Raobs = Raobs(i,j,k);

fit_eqn = fitValues.fit_SS_eqn_sprintf;
fit_eqn = sprintf(fit_eqn, repmat(Raobs, fitValues.numTerms,1));

% Use matrix division: 
X = zeros(length(msat), fitValues.numTerms); 

for i = 1:fitValues.numTerms
    X(:, i) = B1_ref.^i;
end 

try
    fitvals = X \ msat; 
    M0b = fitvals(1);
catch
    disp('An error occurred during matrix division:');
    disp('B1_ref:');
    disp(B1_ref);
    disp('msat values:');
    disp(msat);
    disp('Fit Equation (post sprintf):');
    disp(fit_eqn);
    disp('Matrix X:');
    disp(X);
    return;
end








% Previous version 
opts = fitoptions( 'Method', 'NonlinearLeastSquares','Upper',0.5,'Lower',0.0,'StartPoint',0.1);
opts.Robust = 'Bisquare';

myfittype = fittype( fit_eqn ,'dependent', {'z'}, 'independent',{'b1'},'coefficients', {'M0b'});

disp('B1_ref:');
disp(B1_ref);
disp('Size B1_ref: ');
disp(size(B1_ref));
disp('msat values:');
disp(msat);
disp('size msat: ');
disp(size(msat));
disp('Fit Equation (post sprintf):');
disp(fit_eqn);
disp('Fit Type:');
disp(myfittype);
disp('Fit Options:');
disp(opts);

try
fitpos = fit(B1_ref, msat, myfittype,opts);
catch ME
    % Display the error message
    disp('An error occurred during fitting:');
    disp(ME.message);
    
    % Optional: Display additional details about where the error occurred
    disp('Error identifier:');
    disp(ME.identifier);
    
    disp('Error stack trace:');
    for k = 1:length(ME.stack)
        disp(['In file: ', ME.stack(k).file]);
        disp(['Function: ', ME.stack(k).name]);
        disp(['Line: ', num2str(ME.stack(k).line)]);
    end
end





