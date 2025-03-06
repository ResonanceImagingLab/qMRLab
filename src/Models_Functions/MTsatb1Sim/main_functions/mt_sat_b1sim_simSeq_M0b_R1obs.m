%% Simulate sequence and generate fitequation to cover the spectrum of MTsat
% results for varying B1rms, R1obs and M0b. 
% Please consult README document first to be sure you have downloaded all
% necessary packages. 


function mt_sat_b1sim_simSeq_M0b_R1obs(obj)


% OutputDir = directory where results will be saved 
OutputDir = obj.options.SequenceSimulations_OutputDirectory;

%turboF = [8,80,200];
b1 = linspace(0, 1.5, 30); % relative B1 field
M0b = 0:0.03:0.18; 
T1obs = horzcat(0.6:0.075:2,2.1:0.4:3); %600ms to 4500ms to cover WM to CSF. 
Raobs = 1./T1obs;

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

% Took about 2.5 hours to run 
    tic
    %clear Params outputSamplingTable;

    Params.M0a = obj.Prot.TissueParams.Mat(1);
    Params.Raobs = obj.Prot.TissueParams.Mat(2);
    Params.R = obj.Prot.TissueParams.Mat(3);
    Params.T2a = obj.Prot.TissueParams.Mat(4)/1000;
    Params.R1b = obj.Prot.TissueParams.Mat(5);
    Params.T2b = obj.Prot.TissueParams.Mat(6)/1e6;
    Params.M0b = obj.Prot.TissueParams.Mat(7);

    Params.MTC = obj.Prot.PulseSequenceParams.Mat(1);
    Params.delta = obj.Prot.PulseSequenceParams.Mat(2);  
    Params.flipAngle = obj.Prot.PulseSequenceParams.Mat(3);
    Params.TR = obj.Prot.PulseSequenceParams.Mat(4)/1000;
    Params.numSatPulse = obj.Prot.PulseSequenceParams.Mat(5);
    Params.pulseDur = obj.Prot.PulseSequenceParams.Mat(6)/1000;
    Params.satFlipAngle = obj.Prot.PulseSequenceParams.Mat(7);

    Params.SatPulseShape = obj.options.SequenceSimulations_SatPulseShape;

    Params = mt_sat_b1sim_getSeqParams(Params);

    % Loop variables:
    Params.M0b =  []; 
    Params.Raobs = [];
    Params.Ra = [];

    GRE_sigs = zeros(size(b1,2),size(M0b,2),size(Raobs,2));

    tic
    disp('Simulation starting, may take 2-3 hours to run')
    for i = 1:size(b1,2) % took nearly 5 hours for matrix 25x41x33.

        if (wait)
            % Update waitbar
            ww = ww+1;
            waitbar(ww/size(b1,2), h, sprintf('Simulating with B1+ = %.2f', b1(i)));
        end
    
        for j = 1:size(M0b,2)
            Params.M0b = M0b(j);
            
            for k = 1:size(Raobs,2)
                Params.Raobs = Raobs(k);
                
                GRE_sigs(i,j,k) = mt_sat_b1sim_blochSimFlashSequence_v2(Params,'freqPattern','single', 'satFlipAngle', b1(i)*Params.satFlipAngle);       

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
    
    disp('Simulations Finished. Total time: ')
    toc
    disp('Now Calculating Values...')

    % Delete waitbar
    delete(h);

    %% MTsat calculation
    %reformat Aapp and R1app matrices for 3D calculation
    Aapp = ones(size(GRE_sigs));
    T1app = repmat(T1obs,[7,1,size(b1,2)]);
    T1app = permute(T1app,[3,1,2]);
    

    % Repurpose code from ihMT module
    MTsat_sim_S = ihMT_calcMTsatThruLookupTablewithDummyV3( GRE_sigs,   [], T1app* 1000, [], Aapp,...
        0, Params.numExcitation, Params.TR * 1000, Params.flipAngle, 0);

    MTsatValue_fn = fullfile(OutputDir, strcat('MTsat_sim_S','.mat')); 
    save(MTsatValue_fn,'MTsat_sim_S') % MTsat_sim_S, out{}m '.mat'

 
    %% Clean up then fit:
    MTsat_sim_S(MTsat_sim_S < 0) = NaN;
 
    % Single % Repurpose code from ihMT module
    [fit_SS_eqn, fit_SS_eqn_sprintf, fit_SSsat, numTerms] = ihMT_generateFitequationsV2(M0b(1:6), b1, Raobs, MTsat_sim_S(:,1:6,:));

    % put into one variable for export
    fitValues.fitvals_coeff = fit_SSsat.Coefficients;
    fitValues.fit_SS_eqn = fit_SS_eqn;
    fitValues.fit_SS_eqn_sprintf = fit_SS_eqn_sprintf;
    fitValues.Params = Params; % export params to reference later if desired
    fitValues.numTerms = numTerms; % for fitting later...

    fitValue_fn = fullfile(OutputDir, strcat('fitValues_S','.mat')); 
    save(fitValue_fn,'fitValues')

    img_fn = fullfile(OutputDir, strcat('simFig_S','.png')); 
    ihMT_generateFitSimFigures(M0b(1:6), b1, Raobs, MTsat_sim_S(:,1:6,:), fit_SS_eqn, img_fn)



    toc
%end







