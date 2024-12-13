classdef ihMT < AbstractModel

    % ihMT:   inhomogenuous Magnetization Transfer
    %
    % Assumptions: 
    %         B1+ corrected MT saturation maps taking into account for the
    %         B1+ inhomogeneities effects on the excitation and saturation
    %         pulses.
    %
    % Inputs:
    %   dual               MT-weigthed data. Dual frequency preparation
    %                      pulse.
    %   pos                MT-weigthed data. Positive single sided frequency 
    %                      preparation pulse.
    %   neg                MT-weighted dats. Negative single sided frequency
    %                      preparation pulses 
    %   T1map              T1-weighted data.
    %   M0map              PD-weighted data.
    %   B1                 Normalized transmit excitation field map (B1+). B1+ is defined 
    %                      as a  normalized multiplicative factor such that:
    %                      FA_actual = B1+ * FA_nominal.
    %   (mask)             Binary mask (Accelerate fitting - Optional).
    %            
    % Outputs:
    %
    % Protocol:	
    %
    % Options:
    %   See:
    %       Model.options (general options)
    %       Model.options.Sequencesimulation (to change parameters of the sequence simulation)
    %   
    % Example of command line usage
    %
    %   For more examples: <a href="matlab: qMRusage(ihMT);">qMRusage(ihMT)</a>
    %
    % Author:
    %   Christopher D. Rowley, 2023 (@christopherrowley, @TardifLab - GitHub)
    % Adapted to qMRLab by:
    %   Amie Demmans, 2024
    %
    % References:
    %   Please cite the following if you use this module:
    %     Rowley C.D., Campbell J.S.W., Wu Z., Leppert I.R., Nelson M.C.,
    %     Pike G.B., Tardif C.L. (2023), Optimization of acquisition parameters for 
    %     cortical inhomogeneous magnetization transfer (ihMT) imaging using a rapid 
    %     gradient echo readout. Magn Reson
    %     Med 90(5):1762-1775. doi:10.1002/mrm.29754
    %   In addition to citing the package:
    %     Karakuzu A., Boudreau M., Duval T.,Boshkovski T., Leppert I.R., Cabana J.F., 
    %     Gagnon I., Beliveau P., Pike G.B., Cohen-Adad J., Stikov N. (2020), qMRLab: 
    %     Quantitative MRI analysis, under one umbrella doi: 10.21105/joss.02343

    properties (Hidden=true)
        %onlineData_url = 'https://osf.io/3s9xe/download?version=2';
        % Need to figure out what exact data needs to go here/ which images
        % Need to add atlas images--> create issure on qMRLab 
        % https://github.com/qMRLab/qMRLab/wiki/Guideline:-Uploading-sample-data
    end

    properties 
        MRIinputs = {'dual', 'pos', 'neg', 'T1map', 'M0map', 'b1', 'mask'};
        xnames = {};
        voxelwise = 0; % 0, if the analysis is done matricially
        % 1, if the analysis is done voxel per voxel

        % PulseSequenceParams & Tissue Params
        Prot = struct('PulseSequenceParams', struct('Format',{{'MTC'; 'delta' ; 'flipAngle' ; 'TR(ms)' ; 'numSatPulse' ; 'TurboFactor' ; 'pulseDur(ms)' ; 'satFlipAngle' ; ...
                                     'pulseGapDur(ms)' ; 'DummyEcho' ; 'LowDutyCycle'; 'satTrainPerBoost' ; 'TR_MT(ms)'; 'echoSpacing(ms)'}}, ...
                                    'Mat', [1; 8000; 7; 1140; 6; 80; 0.768; 136; 0.3; 2; 1; 9; 60; 7.66]), ...
                      'TissueParams', struct('Format',{{'M0a'; 'Raobs'; 'R'; 'T2a(ms)'; 'T1D(ms)'; 'R1b'; 'T2b(μs)'; 'M0b'; 'D'}}, ...
                      'Mat', [ 1; 1/1.4; 50; 50; 0.75; 0.25; 11.5; 0.071; 0.8e-3/1e6]));

        % Params.numExcitation is based on TurboFactor and Dummy Echo 
        % freqpattern ill make into dropdown in buttons probably b/c I
        % can't put it in Prot unless it is a number
        % CR_getSeqParams_3prot.m

        % Sequence Simulation fitVals
        fitValues_dual = [];
        fitValues_single = [];

        % R1vsM0b correlation fitVals
        fitValues_Dual = [];
        fitValues_SP = [];
        fitValues_SN = [];

        buttons = {'PANEL', 'SequenceSimulations',7,...
            'AtlasDirectory', 0 ,... 
            'OutputDirectory',0,...
            'TissueType',{'GM','WM'}, ...
            'B0', {'3', '7', '1.5'}, ...
            'FreqPattern',{'dualAlternate','dualContinuous'},...
            'SatPulseShape', {'gausshann', 'gaussian', 'fermi'},...
            'Run Sequence Simulations','pushbutton',...
            'PANEL', 'R1vsM0b Mapping',3,...
            'SeqSimDirectory', 0,...
            'RunR1vsM0bCorrelation',true,...
            'Select Appropriate Directories','pushbutton'};

        options = struct(); 
        previousOptions = struct()

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
                ~isequal(obj.options.SequenceSimulations_AtlasDirectory, obj.previousOptions.SequenceSimulations_AtlasDirectory)||...
                ~isequal(obj.options.SequenceSimulations_OutputDirectory, obj.previousOptions.SequenceSimulations_OutputDirectory)||...
                ~isequal(obj.options.SequenceSimulations_FreqPattern, obj.previousOptions.SequenceSimulations_FreqPattern)||...
                ~isequal(obj.options.SequenceSimulations_SatPulseShape, obj.previousOptions.SequenceSimulations_SatPulseShape)||...
                ~isequal(obj.options.R1vsM0bMapping_SeqSimDirectory, obj.previousOptions.R1vsM0bMapping_SeqSimDirectory)||...
                ~isequal(obj.options.R1vsM0bMapping_RunR1vsM0bCorrelation, obj.previousOptions.R1vsM0bMapping_RunR1vsM0bCorrelation))   
            checkfields = 1; 
        elseif (~isequal(obj.options.SequenceSimulations_RunSequenceSimulations, obj.previousOptions.SequenceSimulations_RunSequenceSimulations)||...
                 ~isequal(obj.options.R1vsM0bMapping_SelectAppropriateDirectories, obj.previousOptions.R1vsM0bMapping_SelectAppropriateDirectories))
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
            obj.Prot.TissueParams.Mat = [Params.M0a, Params.Raobs, Params.R, Params.T2a*1000, ...
                                        Params.T1D*1000, Params.R1b, Params.T2b*1e6, ...
                                        Params.M0b, Params.D]';

        elseif obj.checkupdatedfields == 2
         
            if obj.options.SequenceSimulations_RunSequenceSimulations
                disp('Select directory with atlas images')
                obj.options.SequenceSimulations_AtlasDirectory = uigetdir(pwd);
                disp('Select directory to save fit values')
                obj.options.SequenceSimulations_OutputDirectory = uigetdir(pwd); 
                
                ihMT_simSeq_M0b_R1obs_3prot(obj); 
    
            elseif obj.options.R1vsM0bMapping_SelectAppropriateDirectories
                disp('Select directory where you want values saved')
                obj.options.R1vsM0bMapping_SeqSimDirectory = uigetdir(pwd);
               
                if ~obj.options.R1vsM0bMapping_RunR1vsM0bCorrelation
                    disp('Load fitValues_Dual.mat')
                    [FileName_Dual, PathName_Dual] = uigetfile('*.mat');
                    disp('Load fitValues_SP.mat'); 
                    [FileName_SP, PathName_SP] = uigetfile('*.mat'); 
                    disp('Load fitValues_SN.mat');
                    [FileName_SN, PathName_SN] = uigetfile('*.mat'); 

                    obj.fitValues_Dual = load([PathName_Dual filesep FileName_Dual]); 
                    obj.fitValues_SP = load([PathName_SP filesep FileName_SP]); 
                    obj.fitValues_SN = load([PathName_SN filesep FileName_SN]); 
                else
                    disp('Load dual fit values')
                    [FileName_dual,PathName_dual] = uigetfile('*.mat');
                    disp('Load single fit values')
                    [FileName_single,PathName_single] = uigetfile('*.mat');
                         
                    obj.fitValues_dual = load([PathName_dual filesep FileName_dual]);
                    obj.fitValues_single = load([PathName_single filesep FileName_single]); 

                end
                
            end 

        end 

    end 

    function FitResult = fit(obj,data)
         
         flipA = obj.Prot.PulseSequenceParams.Mat(3);
         TR = obj.Prot.PulseSequenceParams.Mat(4); % ms
         DummyEcho = obj.Prot.PulseSequenceParams.Mat(10);
         echoSpacing = obj.Prot.PulseSequenceParams.Mat(14); % ms 
         numExcitation = obj.Prot.PulseSequenceParams.Mat(6) + DummyEcho;

        if obj.options.R1vsM0bMapping_RunR1vsM0bCorrelation % If box is checked, run correlation 
            fitValues_dual = obj.fitValues_dual;
            fitValues_single = obj.fitValues_single; 
           
            [fitValues_Dual, fitValues_SP, fitValues_SN] = ihMT_R1vsM0b_correlation(obj, data, fitValues_dual, fitValues_single, flipA, TR, DummyEcho, echoSpacing, numExcitation);
        
        else
            fitValues_Dual = obj.fitValues_Dual;
            fitValues_SP = obj.fitValues_SP;
            fitValues_SN = obj.fitValues_SN;
        end

        [sat_dual_c, sat_pos_c, sat_neg_c, ihmt_c] = ihMT_correctMTsat_3prot(data, fitValues_Dual, fitValues_SP, fitValues_SN, flipA, TR, DummyEcho, echoSpacing, numExcitation);
        FitResult.sat_dual_c = sat_dual_c;
        FitResult.sat_pos_c = sat_pos_c; 
        FitResult.sat_neg_c = sat_neg_c; 
        FitResult.ihmt_c = ihmt_c; 
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
















