classdef ihMT < AbstractModel

    properties (Hidden=true)
        onlineData_url = 'https://osf.io/3s9xe/download?version=2';
        % Need to figure out what exact data needs to go here/ which images
        % Need to add atlas images--> create issure on qMRLab 
        % https://github.com/qMRLab/qMRLab/wiki/Guideline:-Uploading-sample-data
    end

    properties 
        MRIinputs = {'MTw_dual', 'MTw_single', 'T1w', 'PDw', 'B1map', 'Mask'};
        xnames = {};
        voxelwise = 0; % 0, if the analysis is done matricially
        % 1, if the analysis is done voxel per voxel
        % Protocol
        Prot = struct('MtParams', struct('Format',{{'delta' ; 'flipAngle' ; 'TR' ; 'numSatPulse' ; 'turboFactor' ; 'pulseDur' ; 'b1' ; ...
                                    'satTrain/Boost' ; 'TR_MT' ; 'pulseGapDur' ; 'DummyEcho' ; 'boosted'}}, ...
                                    'Mat', [8000; 7; 1.14; 6; 80; 0.000768; 11.6; 9; 0.06; 0.0003; 2; 1]), ...
                      'DefaultTissueParams', struct('Format',{{'M0a'; 'Raobs'; 'R'; 'T2a'; 'T1D'; 'R1b'; 'T2b'; 'Ra'; 'M0b'; 'D'}}, ...
                      'Mat', [1; 1/1.4; 50; 50e-3; 7.5e-4; 0.25; 11.5e-6; 1; 0.071; 0.8e-3/1e6]));

        % Params.numExcitation is based on TurboFactor and Dummy Echo 
        % freqpattern ill make into dropdown in buttons probably b/c I
        % can't put it in Prot unless it is a number
        % CR_getSeqParams_3prot.m


        % 
        % ProtStyle = struct('prot_namespace',{{'MTw_dual', 'MTw_single', 'T1w','PDw'}}, ...
        % 'style',repmat({'TableNoButton'},[1,4]));
    
        % fitValues_dual = load('fitValues_dualAlt.mat');
        % fitValues_single = load('fitValues_single.mat');

        buttons = {'DataDirectory', 0 ,... 
            ' OutputDirectory',0,...
            'TissueType',{'WM','GM'}, ...
            'B0', {'3', '7', '1.5'}, ...
            'FreqPattern',{'single','dualAlternate','dualContinuous'},...
            'Run Sequence Simulations','pushbutton'};
       

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

        options = struct(); 

    end 

methods
    function obj = ihMT()
        obj.options = button2opts(obj.buttons);
    end 
end

end 
























