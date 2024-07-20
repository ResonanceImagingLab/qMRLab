%% Simulate sequence and generate fitequation to cover the spectrum of MTsat
% results for varying B1rms, R1obs and M0b. 
% Please consult README document first to be sure you have downloaded all
% necessary packages. 


function ihMT_simSeq_M0b_R1obs_3prot(obj)

%DATADIR = 'Image/Directory';
%DATADIR = '/Users/amiedemmans/Desktop/GitHubOptimizeIHMTimaging/kspaceWeighting/Atlas_reference/';

%load( strcat( '/Users/amiedemmans/Desktop/GitHub/OptimizeIHMTimaging/kspaceWeighting/Atlas_reference/','GM_seg_MNI_152_image.mat'))
load(strcat(obj.options.SequenceSimulations_DataDirectory, '/GM_seg_MNI_152_image.mat'))
load(strcat(obj.options.SequenceSimulations_DataDirectory,'/GM_seg_MNI_152_kspace.mat'))

%OutputDir =  'directory/outputs';
%OutputDir = '/Users/amiedemmans/Desktop/GitHub/OptimizeIHMTimaging/b1Correction/SeqSim2';
OutputDir = obj.options.SequenceSimulations_OutputDirectory;

%turboF = [8,80,200];
b1 = 0:0.75:18;
M0b = 0:0.03:0.18; 
T1obs = horzcat(0.6:0.075:2,2.1:0.4:3); %600ms to 4500ms to cover WM to CSF. 
Raobs = 1./T1obs;

%for z = 1:length(turboF)

% Took about 2.5 hours to run 
    tic
    %clear Params outputSamplingTable;

    Params.M0a = obj.Prot.TissueParams.Mat(1);
    Params.Raobs = obj.Prot.TissueParams.Mat(2);
    Params.R = obj.Prot.TissueParams.Mat(3);
    Params.T2a = obj.Prot.TissueParams.Mat(4);
    Params.T1D = obj.Prot.TissueParams.Mat(5);
    Params.R1b = obj.Prot.TissueParams.Mat(6);
    Params.T2b = obj.Prot.TissueParams.Mat(7);
    Params.M0b = obj.Prot.TissueParams.Mat(8);
    Params.D = obj.Prot.TissueParams.Mat(9);

    Params.MTC = obj.Prot.PulseSequenceParams.Mat(1);
    Params.delta = obj.Prot.PulseSequenceParams.Mat(2);  % delta 
    Params.flipAngle = obj.Prot.PulseSequenceParams.Mat(3);
    Params.TR = obj.Prot.PulseSequenceParams.Mat(4);
    Params.numSatPulse = obj.Prot.PulseSequenceParams.Mat(5);
    Params.TurboFactor = obj.Prot.PulseSequenceParams.Mat(6);
    Params.pulseDur = obj.Prot.PulseSequenceParams.Mat(7);
    Params.satFlipAngle = obj.Prot.PulseSequenceParams.Mat(8);
    Params.pulseGapDur = obj.Prot.PulseSequenceParams.Mat(9);
    Params.DummyEcho = obj.Prot.PulseSequenceParams.Mat(10);
    Params.boosted = obj.Prot.PulseSequenceParams.Mat(11);
    Params.satTrainPerBoost = obj.Prot.PulseSequenceParams.Mat(12);
    Params.TR_MT = obj.Prot.PulseSequenceParams.Mat(13);

    [Params, outputSamplingTable] = ihMT_getSeqParams_3prot(Params);


    % Loop variables:
    Params.M0b =  []; % going to loop over this
    Params.Raobs = [];
    Params.Ra = [];
    %Params.satFlipAngle = Params.pulseDur*360*42.58; 

     gm_m = brain_m;
     fft_gm_m = fft_brain_m;

    GRE_sigd = zeros(size(b1,2),size(M0b,2),size(Raobs,2));
    GRE_sigs = zeros(size(b1,2),size(M0b,2),size(Raobs,2));

    tic
    for i = 1:size(b1,2) % took nearly 5 hours for matrix 25x41x33.
        %Params.b1 = b1(i);
    
        for j = 1:size(M0b,2)
            Params.M0b = M0b(j);
            
            for k = 1:size(Raobs,2)
                Params.Raobs = Raobs(k);
                temp = ihMT_blochSimFlashSequence_v2(Params,'freqPattern','single', 'satFlipAngle', b1(i)*Params.satFlipAngle);       
                GRE_sigs(i,j,k) = ihMT_generate_BSF_scaling_v1( temp, Params, outputSamplingTable, gm_m, fft_gm_m) ;

                if strcmp(obj.options.SequenceSimulations_FreqPattern, 'dualAlternate')
                    temp = ihMT_blochSimFlashSequence_v2(Params,'freqPattern','dualAlternate', 'satFlipAngle', b1(i)*Params.satFlipAngle);
                elseif strcmp(obj.options.SequenceSimulations_FreqPattern, 'dualContinuous')
                    temp = ihMT_blochSimFlashSequence_v2(Params,'freqPattern','dualContinuous', 'satFlipAngle', b1(i)*Params.satFlipAngle);
                end

                GRE_sigd(i,j,k) = ihMT_generate_BSF_scaling_v1( temp, Params, outputSamplingTable, gm_m, fft_gm_m) ;
            end
        end
        disp(i/size(b1,2) *100)  % print percent done...
        toc
    end


    %% MTsat calculation
    %reformat Aapp and R1app matrices for 3D calculation
    Aapp = ones(size(GRE_sigd));
    T1app = repmat(T1obs,[7,1,size(b1,2)]);
    T1app = permute(T1app,[3,1,2]);
    
    flip_rad = Params.flipAngle*pi/180 ; % use the nominal value here 

    MTsat_sim_S = ihMT_calcMTsatThruLookupTablewithDummyV3( GRE_sigs,   [], T1app* 1000, [], Aapp,...
        Params.echoSpacing * 1000, Params.numExcitation, Params.TR * 1000, Params.flipAngle, Params.DummyEcho);
    MTsat_sim_D = ihMT_calcMTsatThruLookupTablewithDummyV3( GRE_sigd,   [], T1app* 1000, [], Aapp,...
        Params.echoSpacing * 1000, Params.numExcitation, Params.TR * 1000, Params.flipAngle, Params.DummyEcho);



    MTsatValue_fn = fullfile(OutputDir, strcat('MTsat_sim_S','.mat')); 
    save(MTsatValue_fn,'MTsat_sim_S') % MTsat_sim_S, out{}m '.mat'

    MTsatValue_fn = fullfile(OutputDir, strcat('MTsat_sim_D','.mat')); 
    save(MTsatValue_fn,'MTsat_sim_D')



    %% Clean up then fit:
    MTsat_sim_S(MTsat_sim_S < 0) = NaN;
    MTsat_sim_D(MTsat_sim_D < 0) = NaN;

    % Single
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



    % Dual
    [fit_SS_eqn, fit_SS_eqn_sprintf, fit_SSsat, numTerms] = ihMT_generateFitequationsV2(M0b(1:6), b1, Raobs, MTsat_sim_D(:,1:6,:));

    % put into one variable for export
    fitValues.fitvals_coeff = fit_SSsat.Coefficients;
    fitValues.fit_SS_eqn = fit_SS_eqn;
    fitValues.fit_SS_eqn_sprintf = fit_SS_eqn_sprintf;
    fitValues.Params = Params; % export params to reference later if desired
    fitValues.numTerms = numTerms; % for fitting later...

    fitValue_fn = fullfile(OutputDir, strcat('fitValues_D','.mat'));
    save(fitValue_fn,'fitValues')

    img_fn = fullfile(OutputDir, strcat('simFig_D','.png')); 
    ihMT_generateFitSimFigures( M0b(1:6), b1, Raobs, MTsat_sim_D(:,1:6,:), fit_SS_eqn, img_fn)


    toc
%end







