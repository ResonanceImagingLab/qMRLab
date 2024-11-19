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

% May need to use non sprintf --> look into both
% Evaluating linear text equations with matlab 
% Look into symbolic toolbox documentation 

fit_eqn = fitValues.fit_SS_eqn;
% fit_eqn = sprintf(fit_eqn, repmat(Raobs, fitValues.numTerms,1));

% Initialize degrees
B1_degree = 0;
Raobs_degree = 0;

% Extract powers from fit_eqn
B1_powers = regexp(fit_eqn, 'b1\.\^(\d+)', 'tokens');
Raobs_powers = regexp(fit_eqn, 'Raobs\.\^(\d+)', 'tokens');

% Extract constant from fit_eqn 
constants =  regexp(fit_eqn, '[\+\-]?\d+\.\d+', 'match');
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
% Handle edge cases
% Step through to see where complex value is coming from 
fitV(imag(fitV)~=0)= NaN;

fitV(fitV<0) = NaN;
[~,temp] = min(abs(fitV-0.1));
M0b = fitV(temp);
if isnan(M0b)
    M0b = 0; 
end 


% fit_eqn = fitValues.fit_SS_eqn_sprintf;
% fit_eqn = sprintf(fit_eqn, repmat(Raobs, fitValues.numTerms,1));
% 
% % Construct vandermonde matrix for matrix division: 
% V = zeros(length(B1_ref), fitValues.numTerms); 
% 
% for i = 1:fitValues.numTerms
%     V(:, i) = B1_ref.^(i-1); % A(i,j)=v(i)(N-j)
% 
% end 
% 
% try
%     fitvals = V \ msat; 
%     M0b = fitvals(1);
% 
% catch
%     disp('An error occurred during matrix division:');
%     disp('B1_ref:');
%     disp(B1_ref);
%     disp('msat values:');
%     disp(msat);
%     disp('Raobs');
%     disp(Raobs);
%     disp('Matrix V:');
%     disp(V);
%     return;
% end

% disp('size b1_ref: ')
% disp(size(B1_ref))
% 
% disp('size msat '); 
% disp(size(msat));

% 

% disp('Fit equation from sprintf: ')
% disp(fit_eqn)


% opts = fitoptions( 'Method', 'NonlinearLeastSquares','Upper',0.5,'Lower',0.0,'StartPoint',0.1);
% opts.Robust = 'Bisquare';
% 
% myfittype = fittype( fit_eqn ,'dependent', {'z'}, 'independent',{'b1'},'coefficients', {'M0b'}); 
% 
% 
% fitpos = fit(B1_ref, msat, myfittype,opts); % insert a try+catch to report values. In the catch -> break. 
% try
%     fitpos = fit(B1_ref, msat, myfittype,opts); % insert a try+catch to report values. In the catch -> break. 
% catch 
%    disp('B1_ref:');
%    disp(B1_ref);
%    disp('msat values:');
%    disp(msat);
%    disp('Fit Equation (post sprintf):');
%    disp(fit_eqn);
%    disp('Fit Type:');
%    disp(myfittype);
%    disp('Fit Options:');
%    disp(opts);
%    return;
% end 
% 
% fitvals = coeffvalues(fitpos);
% 
% M0b = fitvals(1);

%% Calculate Residuals
% Uncomment if you would like to calculate residuals 
% Not including it right now to save coding time 
% solve the equation for the 
% comb_resid = 0;
% for i = 1:size(msat,1)
%     b1 = B1_ref(i);
%     tmp = eval(fit_eqn);
%     resid = abs(tmp - msat(i));
%     comb_resid = comb_resid + resid;
% end


 
%% you can plot to check the fit
% b1_ref = 0:0.25:11;
% msat_calc = fitpos(b1_ref);
% figure;
% plot(b1_ref,msat_calc,'LineWidth',2)
% hold on
% scatter(B1_ref, msat,40,'filled')
%     ax = gca;
%     ax.FontSize = 20; 
%     xlabel('B_{1RMS} (\muT) ', 'FontSize', 20, 'FontWeight', 'bold')
%     ylabel('MT_{sat}', 'FontSize', 20, 'FontWeight', 'bold')
%      %   colorbar('off')
%     legend('hide')
%     text(6.2, 0.0015, strcat('M_{0,app}^B = ',num2str(M0b,'%.3g')), 'FontSize', 16); 
%     ylim([-0.001 20e-3])

    
    
   %% 2k code. 
    
%   b1_ref = 0:0.25:5;
% msat_calc = fitpos(b1_ref);
% figure;
% plot(b1_ref,msat_calc,'LineWidth',2)
% hold on
% scatter(B1_ref, msat,40,'filled')
%     ax = gca;
%     ax.FontSize = 20; 
%     xlabel('B_{1RMS} (\muT) ', 'FontSize', 20, 'FontWeight', 'bold')
%     ylabel('MT_{sat}', 'FontSize', 20, 'FontWeight', 'bold')
%      %   colorbar('off')
%     legend('hide')
%     text(2.2, 0.0015, strcat('M_{0,app}^B = ',num2str(M0b,'%.3g')), 'FontSize', 16); 
%     ylim([-0.001 25e-3])  
    
    
%% Something I was trying 
% fit_eqn = fitValues.fit_SS_eqn;
% % fit_eqn = sprintf(fit_eqn, repmat(Raobs, fitValues.numTerms,1));
% 
% % Initialize degrees
% M0b_degree = 0; 
% B1_degree = 0;
% Raobs_degree = 0;
% 
% % Extract powers from fit_eqn
% M0b_powers = regexp(fit_eqn, 'M0b\.\^(\d+)', 'tokens');
% B1_powers = regexp(fit_eqn, 'b1\.\^(\d+)', 'tokens');
% Raobs_powers = regexp(fit_eqn, 'Raobs\.\^(\d+)', 'tokens');
% 
% if ~isempty(M0b_powers) 
%     M0b_degree = max(cellfun(@(x) str2double(x), [M0b_powers{:}]));
% end
% if ~isempty(B1_powers)
%     B1_degree = max(cellfun(@(x) str2double(x), [B1_powers{:}]));
% end
% if ~isempty(Raobs_powers)
%     Raobs_degree = max(cellfun(@(x) str2double(x), [Raobs_powers{:}]));
% end 
% 
% % Construct vandermonde matrix for matrix division: 
% V = zeros(length(B1_ref), fitValues.numTerms); 
% 
% % numTerms = possible combinations of powers of M0b, B1 and R1
% idx = 1;
% for i = 0:M0b_degree
%     for j = 0:B1_degree
%         for k = 0:Raobs_degree
%             % The terms of the model will correspond to powers of M0b, b1, and Raobs
%             V(:, idx) = (msat.^(i)) .* (B1_ref.^(j)) .* (Raobs.^(k));
%             idx = idx + 1;
%         end
%     end
% end 
% 
% try
%     fitvals = V \ msat; 
%     M0b = fitvals(1);
% catch
%     disp('An error occurred during matrix division:');
%     disp('B1_ref:');
%     disp(B1_ref);
%     disp('msat values:');
%     disp(msat);
%     disp('Raobs');
%     disp(Raobs);
%     disp('Matrix V:');
%     disp(V);
%     return;
% end    