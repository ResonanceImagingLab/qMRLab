function [M0b, maskval]= ihMT_fit_M0b_v2(B1_ref,Raobs, msat,fitValues)
% V2 corrects for flexibility in the model terms used to fit the simulation
% space. 

% we have 4 data points for saturation level for dual saturation
% each point is collected at a set B1 value and achieves a saturation level

% B1_ref -> contains the x data points. Take the relative B1map and multiply by the nominal B1 values
% Raobs -> is the measured R1 map, from VFA 
% msat -> is the Y data points to fit. MTsat values for each of the X points. 
% fitValues -> is derived from "simWith_M0b_R1obs.m" script

% Exports:
% M0b is the fit value for M0b
% maskval is a mask to signify poor fit regions
% comb_resid == the summation of the residuals from all fit points. 

% %Debugging test parameters
% i = 64; j = 66; k = 46; % Genu of Corpus Callosum for paper.
% % fitValues = fitValues_dual;
% % B1_ref = squeeze(b1_comb_scaled(i,j,k,:)) % this should give us 4 data points
% % msat = squeeze(dual_s(i,j,k,:)) % test values for "z"
% Raobs = R1_s(i,j,k)
% 
% fitValues = fitValues2k;
% B1_ref = squeeze(b1_3p26(i,j,k,:)) % this should give us 4 data points
% msat = msat_irl_2k(i,j,k) % test values for "z"

% Don't need comb_resid for this 
% Re do fit_eqn with anonymous functions??

if ~isfield(fitValues,'numTerms')
    fitValues.numTerms = 90;
end

maskval = 0;

if (min(msat) == 0) || (max(isnan(msat)) > 0) % fit will be poor
    M0b = 0; 
    maskval = 1; % just generate a mask to note regions which might have a bad fit. 
    return;
end

fit_eqn = fitValues.fit_SS_eqn;
% fit_eqn = sprintf(fit_eqn, repmat(Raobs, fitValues.numTerms,1));

% Initialize degrees
B1_degree = 0;
Raobs_degree = 0;

% Extract powers from fit_eqn
B1_powers = regexp(fit_eqn, 'b1\.\^(\d+)', 'tokens');
Raobs_powers = regexp(fit_eqn, 'Raobs\.\^(\d+)', 'tokens');

% Extract constant from fit_eqn 
%tokens =  regexp(fit_eqn, '[\+\-]?\d+\.\d+', 'match');
% modify to handle scientific notation
constants = regexp(fit_eqn, '([+-]?\d*\.?\d+(?:[eE][-+]?\d+)?)(?=\*)', 'match');
constants = str2double(constants); 


if ~isempty(B1_powers)
    B1_degree = (cellfun(@(x) str2double(x), [B1_powers{:}]));
end
if ~isempty(Raobs_powers)
    Raobs_degree = (cellfun(@(x) str2double(x), [Raobs_powers{:}]));
end 

V = zeros(1, 3); 
for j = 1:90
    
        % The terms of the model will correspond to powers of B1, and Raobs
        Value =  constants(j) * (B1_ref.^(B1_degree(j))) .* (Raobs.^(Raobs_degree(j)));
        if j<31
            V(3) = V(3)+Value;
        elseif j>=31 && j<61
            V(2) = V(2)+Value;
        else 
            V(1) = V(1)+Value;
        end 
end

V(3) = V(3)-msat;

fitV = roots(V);

fitV(imag(fitV)~=0)= NaN;
fitV(fitV<0) = NaN;

[~,temp] = min(abs(fitV-0.1));
M0b = fitV(temp);

if isnan(M0b)
    M0b = 0; 
end 


