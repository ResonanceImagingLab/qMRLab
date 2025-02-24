function [T1corr, R1corr, S0corr] = vfa_t1_spoil_calcMaps(data, a1, a2, TR1, TR2,...
    smallFlipApprox, simCoeffs)


% This function ports hMRI functionality to qMRLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(data, 'Mask') || isempty(data.Mask)
    data.Mask = zeros(size(data.LowFlipAngle));
    data.Mask(data.LowFlipAngle>100) =1;
end
        
        %disp(size(data.LowFlipAngle));

% Used to calculate Maps

R1 = vfa_t1_spoil_hmri_calc_R1(struct( 'data', data.LowFlipAngle, 'fa', a1, 'TR', TR1, 'B1', data.B1map),...
                   struct( 'data', data.HighFlipAngle, 'fa', a2, 'TR', TR2, 'B1', data.B1map), smallFlipApprox);

R1 = R1.*data.Mask;

T1 = 1./R1.*data.Mask ; % convert to milliseconds

T1(T1<0) = 0;
T1(isnan(T1)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply spoiling correction:

% Note in Preibisch and Deichmann 2009, A and B are quadratic functions
% dependent on b1map.

% Make correction factors:
A = simCoeffs(1,1)*data.B1map.^2 + simCoeffs(1,2)*data.B1map + simCoeffs(1,3);
B = simCoeffs(2,1)*data.B1map.^2 + simCoeffs(2,2)*data.B1map + simCoeffs(2,3);

% Apply to the T1 map - in milliseconds as simulations are done in ms
T1corr = A + B.*T1;
T1corr = T1corr.*data.Mask;
T1corr(T1corr<0) = 0;
T1corr(isnan(T1corr)) = 0;

R1corr = 1./T1corr;
R1corr = R1corr.*data.Mask;
R1corr(T1corr<0) = 0;
R1corr(isnan(R1corr)) = 0;


% Now recalculate M0;
flip_a = (a1*data.B1map); % correct for B1
x = cos(flip_a) ;
y = exp(-TR1./T1corr);

% Solve for M0 using equation for flass image.
M0_1 = data.LowFlipAngle.*(1-x.*y)./ ( (1-y) .*sin(flip_a));

% second image
flip_a = (a2*data.B1map); % correct for B1
x = cos(flip_a) ;
y = exp(-TR2./T1corr);

% Solve for M0 using equation for flass image.
M0_2 = data.HighFlipAngle.*(1-x.*y)./ ( (1-y) .*sin(flip_a));

S0corr = (M0_1 + M0_2)./2 .* data.Mask; 
S0corr(S0corr<0) = 0;
S0corr(isnan(S0corr)) = 0;
