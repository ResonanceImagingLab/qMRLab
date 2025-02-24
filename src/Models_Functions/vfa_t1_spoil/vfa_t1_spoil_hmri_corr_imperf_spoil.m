function vfa_t1_spoil_hmri_corr_imperf_spoil(obj, saveDir)
%==========================================================================
% PURPOSE
% Compute coefficients to correct for the effect of imperfect spoiling on
% T1 estimation as described in (Preibisch & Deichmann, MRM 2009).
%
% Numerical simulations are performed using the EPG formalism described in
% Malik et al., MRM 2017 and available here:
% https://github.com/mriphysics/EPG-X
%
% The parameters used for the simulation and the resulting correction
% factors are written to a JSON file in the output folder specified by the
% user.
%==========================================================================

disp('--- Calculating Imperfect Spoiling Correction Coefficients ---');
%% ***********************************************%%
% 1./ Numerical simulations with EPG
%*************************************************%%
%%
% Get sequence parameters
FA      =   [obj.Prot.VFAProtocol.Mat(1,1), obj.Prot.VFAProtocol.Mat(2,1)]; % Flip angles [deg]
TR      =   [obj.Prot.VFAProtocol.Mat(1,2), obj.Prot.VFAProtocol.Mat(2,2)]; % [ms]
Phi0    =   obj.Prot.SpoilingCorrection.Mat(1);                                   % [deg]
B1range =   obj.Prot.SpoilingCorrection.Mat(7) :0.05: obj.Prot.SpoilingCorrection.Mat(8); % convert such that 100% = 1
Gdur    =   obj.Prot.SpoilingCorrection.Mat(2);                                   % [ms]
Gamp    =   obj.Prot.SpoilingCorrection.Mat(3);                                   % [mT/m]

if length(Gdur)~= length (Gamp)
    error('The vectors of gradient durations and amplitudes must have the same length!')
end

%% Get tissue parameters
T1range     = obj.Prot.SpoilingCorrection.Mat(4) :100: obj.Prot.SpoilingCorrection.Mat(5);     %[ms]
T2range     = obj.Prot.SpoilingCorrection.Mat(6);     % [ms]
D           = obj.Prot.SpoilingCorrection.Mat(9);   % [um^2/ms]

%% Build structure "diff" to account for diffusion effect
diff        = struct;
diff.D      = D*1e-9;
diff.G      = Gamp;
diff.tau    = Gdur;

%% Run EPG simulation
nT1 = length(T1range);
nT2 = length(T2range);
nB1 = length(B1range);
S1  = zeros([nT1 nT2 nB1]);
S2  = zeros([nT1 nT2 nB1]);
disp('-------- Simulating signals ---------');

for T1val = 1 : nT1 % loop over T1 values, can use parfor for speed
    
    T1 = T1range(T1val);
    npulse = floor(15*T1/min(TR));   % To ensure steady state signal
    
    for T2val = 1 : nT2
        T2 = T2range(T2val);
        
        for B1val = 1 : nB1  % loop over B1+ values
            B1eff = B1range(B1val);
            
            %%% make train of flip angles and their phases
            phi_train = RF_phase_cycle(npulse,Phi0); % phase of the RF pulses
            alpha_train1 = d2r(FA(1).*B1eff)*ones([1 npulse]); % flip angles of the PDw acquisitions
            alpha_train2 = d2r(FA(2).*B1eff)*ones([1 npulse]); % flip angles of the T1w acquisitions
            
            % Calculate signals via EPG:
            
            %PDw
            F0 = EPG_GRE(alpha_train1, phi_train, TR(1), T1, T2, 'diff', diff);
            S1(T1val,T2val,B1val) = (abs(F0(end)));
            %T1w
            F0 = EPG_GRE(alpha_train2, phi_train, TR(2), T1, T2, 'diff', diff);
            S2(T1val,T2val,B1val) = (abs(F0(end)));
            
        end
    end
    disp( [num2str(T1val/nT1 *100), "% done simulations"]);
end
disp('-------- Simulation done. ---------');


%% ***********************************************%%
% 2./ Fitting T1=A(B1eff)+B(B1eff)*T1app
%*************************************************%%
disp('-------- Determining Coefficients. ---------');
ABcoeff = zeros(2, nB1);
T1app = zeros(nB1, nT1, nT2);
for B1val = 1 : nB1
    
    B1eff = B1range(B1val);
    
    % Calculate T1app, accounting for B1+
    T1app(B1val,:,:) = 1./vfa_t1_spoil_hmri_calc_R1(...
        struct('data',S1(:,:,B1val),'fa',d2r(FA(1)),'TR',TR(1),'B1',B1eff),...
        struct('data',S2(:,:,B1val),'fa',d2r(FA(2)),'TR',TR(2),'B1',B1eff),...
        obj.options.IncompleteSpoiling_smallAngleApprox);
    
    % build matrix X with column of ones and column of T1app
    X = ones([nT1*nT2 2]);
    X(:,2) = T1app(B1val,:);
    ABcoeff(:, B1val) = pinv(X)*repmat(T1range, [1 nT2]).';
    
end

%% *********************************************************%%
% 3./ Fitting A=P(B1eff) and B=P(B1eff) with 2nd degree polynom
%***********************************************************%%
polyCoeffA = polyfit(B1range, ABcoeff(1,:), 2);
polyCoeffB = polyfit(B1range, ABcoeff(2,:), 2);


%% *********************************************************%%
% 4./ Compute RMSE on T1app and T1
%***********************************************************%%
T1app = T1app(:,:);
T1corr = repmat(polyval(polyCoeffA, B1range).',[1 nT1*nT2])+ repmat(polyval(polyCoeffB, B1range).',[1 nT1*nT2]).*T1app;
T1_Corr_Err = (T1corr - repmat(T1range, [nB1 nT2]))./repmat(T1range, [nB1 nT2])*100;
T1_App_Err = (T1app - repmat(T1range, [nB1 nT2]))./repmat(T1range, [nB1 nT2])*100;


RMSE_Corr = sqrt(mean(T1_Corr_Err(:).^2));
RMSE_App = sqrt(mean(T1_App_Err(:).^2));


%% *********************************************************%%
% 5./ Write parameters and correction factors in a json file
% in the selected output directory
%***********************************************************%%
P2_a = round(polyCoeffA,4);
P2_b = round(polyCoeffB,4);
RMSE_App=round(RMSE_App,3);
RMSE_Corr=round(RMSE_Corr,3);
if obj.options.IncompleteSpoiling_smallAngleApprox
    stip= 1;
else
    stip = 0;
end

% Output Matrix:
hMRI_spoiling_sims = zeros(3,3);
hMRI_spoiling_sims(1,:) = P2_a;
hMRI_spoiling_sims(2,:) = P2_b;
hMRI_spoiling_sims(3,:) = [RMSE_App, RMSE_Corr, stip];

output_fn = fullfile(saveDir, strcat('hMRI_spoiling_sims','.mat')); 
    save(output_fn,'hMRI_spoiling_sims') % MTsat_sim_S, out{}m '.mat'


% %% *********************************************************%%
% % 5./ Write parameters and correction factors in a json file
% % in the selected output directory
% %***********************************************************%%
% Results.Input = job;
% Results.Output.P2_a = round(polyCoeffA,4);
% Results.Output.P2_b = round(polyCoeffB,4);
% Results.Output.small_angle_approx = job.small_angle_approx;
% Results.Output.RMSE_percent.T1app=round(RMSE_App,3);
% Results.Output.RMSE_percent.T1corr=round(RMSE_Corr,3);
% 
% Results.ToCopy{1}    =['hmri_def.MPMacq_set.names{NN} = ''' job.prot_name ''';' ];
% Results.ToCopy{end+1}=['hmri_def.MPMacq_set.tags{NN}  = ''' strrep(job.prot_name,' ','') ''';'];
% Results.ToCopy{end+1}=['hmri_def.MPMacq_set.vals{NN}  = [' num2str([TR FA]) '];'];
% Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.tag = ''' strrep(job.prot_name,' ','') ''';' ];
% Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.P2_a = [' num2str(round(polyCoeffA,4)) '];'];
% Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.P2_b = [' num2str(round(polyCoeffB,4)) '];'];
% Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.small_angle_approx = ' mat2str(job.small_angle_approx) ';'];
% Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.enabled = hmri_def.imperfectSpoilCorr.enabled;'];
% 
% results_filename = fullfile(job.outdir,[strrep(job.prot_name,' ',''),'.json']);
% 
% spm_jsonwrite(results_filename{1},Results,struct('indent','\t'));

end