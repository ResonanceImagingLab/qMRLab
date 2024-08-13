classdef ihMT < AbstractModel

    properties (Hidden=true)
        onlineData_url = 'https://osf.io/3s9xe/download?version=2';
        % Need to figure out what exact data needs to go here/ which images
        % Need to add atlas images--> create issure on qMRLab 
        % https://github.com/qMRLab/qMRLab/wiki/Guideline:-Uploading-sample-data
    end

    properties 
        MRIinputs = {'MTw_dual', 'MTw_single_pos', 'MTw_single_neg' , 'T1map', 'S0map', 'B1map', 'Mask'};
        xnames = {};
        voxelwise = 0; % 0, if the analysis is done matricially
        % 1, if the analysis is done voxel per voxel
        % Protocol
        Prot = struct('PulseSequenceParams', struct('Format',{{'MTC'; 'delta' ; 'flipAngle' ; 'TR' ; 'numSatPulse' ; 'TurboFactor' ; 'pulseDur' ; 'satFlipAngle' ; ...
                                     'pulseGapDur' ; 'DummyEcho' ; 'LowDutyCycle'; 'satTrainPerBoost' ; 'TR_MT'}}, ...
                                    'Mat', [1; 8000; 7; 1.14; 6; 80; 0.768/1000; 11.6; 0.3/1000; 2; 1; 9; 0.06]), ...
                      'TissueParams', struct('Format',{{'M0a'; 'Raobs'; 'R'; 'T2a'; 'T1D'; 'R1b'; 'T2b'; 'M0b'; 'D'}}, ...
                      'Mat', [ 1; 1/1.4; 50; 50e-3; 7.5e-4; 0.25; 11.5e-6; 0.071; 0.8e-3/1e6]));

        % Params.numExcitation is based on TurboFactor and Dummy Echo 
        % freqpattern ill make into dropdown in buttons probably b/c I
        % can't put it in Prot unless it is a number
        % CR_getSeqParams_3prot.m


        buttons = {'PANEL', 'SequenceSimulations',6,...
            'DataDirectory', 0 ,... 
            'OutputDirectory',0,...
            'TissueType',{'GM','WM'}, ...
            'B0', {'3', '7', '1.5'}, ...
            'FreqPattern',{'dualAlternate','dualContinuous'},...
            'Run Sequence Simulations','pushbutton',...
            'PANEL', 'R1vsM0b Mapping',3,...
            'DataDirectory', 0,...
            'OutputDirectory',0,...
            'Run R1vsM0b Mapping','pushbutton',...
            'PANEL', 'Calculate ihMTsat',3,...
            'DataDirectory', 0,...
            'OutputDirectory',0,...
            'Run ihMTsat Calculation', 'pushbutton'};

        options = struct(); 
        previousOptions = struct();

    end 

methods

    function obj = ihMT()
        obj.options = button2opts(obj.buttons);
        obj.previousOptions = obj.options;
        %obj = UpdateFields(obj);
    end 

    function checkfields = checkupdatedfields(obj)
        if (~isequal(obj.options.SequenceSimulations_TissueType, obj.previousOptions.SequenceSimulations_TissueType) || ...
                ~isequal(obj.options.SequenceSimulations_B0, obj.previousOptions.SequenceSimulations_B0)||...
                ~isequal(obj.options.SequenceSimulations_DataDirectory, obj.previousOptions.SequenceSimulations_DataDirectory)||...
                ~isequal(obj.options.SequenceSimulations_OutputDirectory, obj.previousOptions.SequenceSimulations_OutputDirectory)||...
                ~isequal(obj.options.SequenceSimulations_FreqPattern, obj.previousOptions.SequenceSimulations_FreqPattern)||...
                ~isequal(obj.options.R1vsM0bMapping_DataDirectory, obj.previousOptions.R1vsM0bMapping_DataDirectory)||...
                ~isequal(obj.options.R1vsM0bMapping_OutputDirectory, obj.previousOptions.R1vsM0bMapping_OutputDirectory)||...
                ~isequal(obj.options.CalculateihMTsat_DataDirectory, obj.previousOptions.CalculateihMTsat_DataDirectory)||...
                ~isequal(obj.options.CalculateihMTsat_OutputDirectory, obj.previousOptions.CalculateihMTsat_OutputDirectory))   

            checkfields = 1; 
        elseif (~isequal(obj.options.SequenceSimulations_RunSequenceSimulations, obj.previousOptions.SequenceSimulations_RunSequenceSimulations)||...
                 ~isequal(obj.options.R1vsM0bMapping_RunR1vsM0bMapping, obj.previousOptions.R1vsM0bMapping_RunR1vsM0bMapping)||...
                 ~isequal(obj.options.CalculateihMTsat_RunihMTsatCalculation, obj.previousOptions.CalculateihMTsat_RunihMTsatCalculation))
            checkfields = 2;

        else
            checkfields = 0;
        end 
    end 

    function obj = UpdateFields(obj)

        if obj.checkupdatedfields == 1

            obj.previousOptions = obj.options;

            Params.B0 = str2double(obj.options.SequenceSimulations_B0);
            Params.TissueType = obj.options.SequenceSimulations_TissueType;
            Params = ihMT_defaultCortexTissueParams(Params);
            obj.Prot.TissueParams.Mat = [Params.M0a, Params.Raobs, Params.R, Params.T2a, ...
                                        Params.T1D, Params.R1b, Params.T2b, ...
                                        Params.M0b, Params.D]';
            PulseOpt = ihMT_pulseSeqParams(obj.options);
            obj.Prot.PulseSequenceParams.Mat = [PulseOpt.MTC, PulseOpt.delta, PulseOpt.flipAngle, PulseOpt.TR, PulseOpt.numSatPulse,...
                                            PulseOpt.TurboFactor, PulseOpt.pulseDur, PulseOpt.satFlipAngle, PulseOpt.pulseGapDur, ...
                                            PulseOpt.DummyEcho, PulseOpt.boosted, PulseOpt.satTrainPerBoost, PulseOpt.TR_MT]';
        elseif obj.checkupdatedfields == 2
         
            if obj.options.SequenceSimulations_RunSequenceSimulations
                obj.options.SequenceSimulations_DataDirectory = uigetdir(pwd, 'Select directory where images are');
                obj.options.SequenceSimulations_OutputDirectory = uigetdir(pwd, 'Select directory where you want values saved'); 
    
                ihMT_simSeq_M0b_R1obs_3prot(obj); 
    
            elseif obj.options.R1vsM0bMapping_RunR1vsM0bMapping
                obj.options.R1vsM0bMapping_DataDirectory = uigetdir(pwd, 'Select directory where fit vals are');
                obj.options.R1vsM0bMapping_OutputDirectory = uigetdir(pwd, 'Select directory where you want values saved');
    
            elseif obj.options.CalculateihMTsat_RunihMTsatCalculation
                obj.options.CalculateihMTsat_DataDirectory = uigetdir(pwd, 'Select directory where R1vsM0b vals are');
                obj.options.CalculateihMTsat_OutputDirectory = uigetdir(pwd, 'Select directory where you want values saved');
            end 

        end 

    end 


end

end 



% 
% ProtStyle = struct('prot_namespace',{{'MTw_dual', 'MTw_single', 'T1w','PDw'}}, ...
% 'style',repmat({'TableNoButton'},[1,4]));

% fitValues_dual = load('fitValues_dualAlt.mat');
% fitValues_single = load('fitValues_single.mat');

%Option panel: ModelBasedB1corrected parameters
    % buttons ={'PANEL','Sequence simulation',15,...
    % 'B1rms',9,...
    % 'Number saturation pulse',2,...
    % 'Pulse duration',0.768,...
    % 'Pulse gap duration',0.6,...
    % 'TR',28,...
    % 'WExcDur',3,...
    % 'Number excitation',1,...
    % 'Frequency pattern',{'dualAlternate','single','dualContinuous'},...
    % 'Delta',7000,...
    % 'FlipAngle',9,...
    % 'Saturation pulse shape',{'hanning','gaussian','square'},...
    % '##fitValues Directory',10,...
    % '##fitValues Name',10,...
    % '##MTsatValues Name',10,...
    % 'Run Sequence Simulation','pushbutton',...
    % 'PANEL','Correlate M0bapp VS R1',2,...
    % 'Same Imaging Protocol',true,...
    % 'b1rms',6.8};




% Something I was trying 
    % function PulseOpt = getPulseSeqParams(obj)
    %     obj = ihMT_getSeqParams_3prot(obj);
    % 
    %     PulseOpt.MTC = obj.Prot.PulseSequenceParams.Mat(1);
    %     PulseOpt.delta = obj.Prot.PulseSequenceParams.Mat(2);
    %     PulseOpt.flipAngle = obj.Prot.PulseSequenceParams.Mat(3);
    %     PulseOpt.TR = obj.Prot.PulseSequenceParams.Mat(4);
    %     PulseOpt.numSatPulse = obj.Prot.PulseSequenceParams.Mat(5);
    %     PulseOpt.TurboFactor = obj.Prot.PulseSequenceParams.Mat(6);
    %     PulseOpt.pulseDur = obj.Prot.PulseSequenceParams.Mat(7);
    %     PulseOpt.satFlipAngle = obj.Prot.PulseSequenceParams.Mat(8);
    %     PulseOpt.pulseGapDur = obj.Prot.PulseSequenceParams.Mat(9);
    %     PulseOpt.DummyEcho = obj.Prot.PulseSequenceParams.Mat(10);
    %     PulseOpt.boosted = obj.Prot.PulseSequenceParams.Mat(11);
    %     PulseOpt.satTrainPerBoost = obj.Prot.PulseSequenceParams.Mat(12);
    %     PulseOpt.TR_MT = obj.Prot.PulseSequenceParams.Mat(13);
    % end
















