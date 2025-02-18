classdef mt_sat_b1sim < AbstractModel
% mt_sat_b1sim: inhomogenuous Magnetization Transfer
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
%   b1                 Normalized transmit excitation field map (B1+).
%   (mask)             Binary mask. 
% Outputs:
%   sat_dual_c         Corrected MT-weighted data for dual frequency
%                      preparation pulse.
%   sat_pos_c          Corrected MT-weighted data for positive single
%                      sided frequency preparation pulse.
%   sat_neg_c          Corrected MT-weighted data for negative single
%                      sided frequency preparation pulse.
%   mt_sat_b1sim_c             B1-corrected mt_sat_b1sim image.
% Protocol:	
%   PulseSequenceParams  Default pulse sequence parameters 
%   TissueParams         Default tissue parameters 
% Options:
%   See:
%       Model.options (general options)
%       Model.options.Sequencesimulation (to change parameters of the sequence simulation)
% Author:
%   Christopher D. Rowley, 2023 (@christopherrowley, @ResonanceImagingLab - GitHub)
% Adapted to qMRLab by:
%   Amie Demmans, 2024 (@amie-demmans, @ResonanceImagingLab - GitHub) 
% References:
%   Please cite the following if you use this module:
%     Rowley, C. D., Campbell, J. S., Wu, Z., Leppert, I. R., Rudko, D. A., Pike, G. B., 
%     & Tardif, C. L. (2021). A model‐based framework for correcting inhomogeneity effects 
%     in magnetization transfer saturation and inhomogeneous magnetization transfer 
%     saturation maps. Magnetic resonance in medicine, 86(4), 2192-2207.4

%   In addition to citing the package:
%     Karakuzu A., Boudreau M., Duval T.,Boshkovski T., Leppert I.R., Cabana J.F., 
%     Gagnon I., Beliveau P., Pike G.B., Cohen-Adad J., Stikov N. (2020), qMRLab: 
%     Quantitative MRI analysis, under one umbrella doi: 10.21105/joss.02343

properties (Hidden=true)
    
end

properties 
    MRIinputs = {'mt-weighted', 'T1map', 'M0map', 'b1', 'mask'};
    xnames = {};
    voxelwise = 0; 

    % PulseSequenceParams & Tissue Params
    Prot = struct('PulseSequenceParams', struct('Format',{{'MTC'; 'delta(Hz)' ; 'flipAngle(deg)' ; 'TR(ms)' ; 'numSatPulse' ; ...
                                 'pulseDur(ms)' ; 'satFlipAngle(deg)' }}, ...
                                'Mat', [1; 2000; 6; 30; 1; 4; 220]), ...
                  'TissueParams', struct('Format',{{'M0a'; 'Raobs'; 'R'; 'T2a(ms)';  'R1b'; 'T2b(μs)'; 'M0b'}}, ...
                  'Mat', [ 1; 1/1.4; 50; 50; 0.25; 11.5; 0.071]));

    % Sequence Simulation fitVals
    fitValues_single = [];

    % R1vsM0b correlation fitVals
    fitValues_SP = [];

    buttons = {'PANEL', 'SequenceSimulations',7,...
        'OutputDirectory',0,...
        'TissueType',{'GM','WM'}, ...
        'B0', {'3', '7', '1.5'}, ...
        'SatPulseShape', {'gaussian', 'gausshann', 'fermi'},...
        'Run Sequence Simulations','pushbutton',...
        'PANEL', 'R1vsM0b Mapping',3,...
        'SeqSimDirectory', 0,...
        'RunR1vsM0bCorrelation',true,...
        'Load fit Value Files','pushbutton'};

    options = struct(); 
    previousOptions = struct()

end 

methods

    function obj = mt_sat_b1sim()
        obj.options = button2opts(obj.buttons);
        obj.previousOptions = obj.options;
        %obj = UpdateFields(obj);
    end 

    function checkfields = checkupdatedfields(obj)
        if (~isequal(obj.options.SequenceSimulations_TissueType, obj.previousOptions.SequenceSimulations_TissueType) || ...
                ~isequal(obj.options.SequenceSimulations_B0, obj.previousOptions.SequenceSimulations_B0)||...
                ~isequal(obj.options.SequenceSimulations_OutputDirectory, obj.previousOptions.SequenceSimulations_OutputDirectory)||...
                ~isequal(obj.options.SequenceSimulations_FreqPattern, obj.previousOptions.SequenceSimulations_FreqPattern)||...
                ~isequal(obj.options.SequenceSimulations_SatPulseShape, obj.previousOptions.SequenceSimulations_SatPulseShape)||...
                ~isequal(obj.options.R1vsM0bMapping_SeqSimDirectory, obj.previousOptions.R1vsM0bMapping_SeqSimDirectory)||...
                ~isequal(obj.options.R1vsM0bMapping_RunR1vsM0bCorrelation, obj.previousOptions.R1vsM0bMapping_RunR1vsM0bCorrelation))   
            checkfields = 1; 
        elseif (~isequal(obj.options.SequenceSimulations_RunSequenceSimulations, obj.previousOptions.SequenceSimulations_RunSequenceSimulations)||...
                 ~isequal(obj.options.R1vsM0bMapping_LoadfitValueFiles, obj.previousOptions.R1vsM0bMapping_LoadfitValueFiles))
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
            Params = mt_sat_b1sim_defaultCortexTissueParams(Params);
            obj.Prot.TissueParams.Mat = [Params.M0a, Params.Raobs, Params.R, Params.T2a*1000, ...
                                        Params.R1b, Params.T2b*1e6, Params.M0b]';

        elseif obj.checkupdatedfields == 2
         
            if obj.options.SequenceSimulations_RunSequenceSimulations
                disp('Select directory to save fit values')
                obj.options.SequenceSimulations_OutputDirectory = uigetdir(pwd); 
                
                mt_sat_b1sim_simSeq_M0b_R1obs(obj); 
    
            elseif obj.options.R1vsM0bMapping_LoadfitValueFiles
 
                if ~obj.options.R1vsM0bMapping_RunR1vsM0bCorrelation
                    disp('Load fitValues_SP.mat'); 
                    [FileName_SP, PathName_SP] = uigetfile('*.mat'); 
                    obj.fitValues_SP = load([PathName_SP filesep FileName_SP]); 
             
                else
                    disp('Select directory where you want values saved')
                    obj.options.R1vsM0bMapping_SeqSimDirectory = uigetdir(pwd);
                    disp('Load simulated MTsat fit values')
                    [FileName_single,PathName_single] = uigetfile('*.mat');
                         
                    obj.fitValues_single = load([PathName_single filesep FileName_single]); 

                end
                
            end 

        end 

    end 

    function FitResult = fit(obj,data)
         
         flipA = obj.Prot.PulseSequenceParams.Mat(3);
         TR = obj.Prot.PulseSequenceParams.Mat(4); % ms
         OutputDir =  obj.options.R1vsM0bMapping_SeqSimDirectory;

        if obj.options.R1vsM0bMapping_RunR1vsM0bCorrelation % If box is checked, run correlation 
            fitValues_single = obj.fitValues_single; 
           
            fitValues_SP = mt_sat_b1sim_R1vsM0b_correlation(data, fitValues_single, flipA, TR, 0, 0, 1, OutputDir);
        
        else
            fitValues_SP = obj.fitValues_SP;     
        end

        sat_pos_c = mt_sat_b1sim_correctMTsat(data, fitValues_SP, flipA, TR, 0, 0, 1);
        FitResult.sat_pos_c = sat_pos_c; 

    end 


end

end 


















