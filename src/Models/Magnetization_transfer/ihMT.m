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
        Prot = struct('MTw_dual',struct('Format',{{'FlipAngle' 'TR'}},...
                                   'Mat',  [9 0.028]),...
                      'MTw_single',struct('Format',{{'FlipAngle' 'TR'}},...
                                   'Mat',  [9 0.028]),...
                      'T1w',struct('Format',{{'FlipAngle' 'TR'}},...
                                   'Mat',  [20 0.030]),...
                      'PDw',struct('Format',{{'FlipAngle' 'TR'}},...
                                   'Mat',  [5 0.030]));
                               
        ProtStyle = struct('prot_namespace',{{'MTw_dual', 'MTw_single', 'T1w','PDw'}}, ...
        'style',repmat({'TableNoButton'},[1,4]));
    
        % fitValues_dual = load('fitValues_dualAlt.mat');
        % fitValues_single = load('fitValues_single.mat');
        
        %Option panel: ModelBasedB1corrected parameters
            buttons ={'PANEL','Sequence simulation',15,...
            'B1rms',9,...
            'Number saturation pulse',2,...
            'Pulse duration',0.768,...
            'Pulse gap duration',0.6,...
            'TR',28,...
            'WExcDur',3,...
            'Number excitation',1,...
            'Frequency pattern',{'dualAlternate','single','dualContinuous'},...
            'Delta',7000,...
            'FlipAngle',9,...
            'Saturation pulse shape',{'hanning','gaussian','square'},...
            '##fitValues Directory',10,...
            '##fitValues Name',10,...
            '##MTsatValues Name',10,...
            'Run Sequence Simulation','pushbutton',...
            'PANEL','Correlate M0bapp VS R1',2,...
            'Same Imaging Protocol',true,...
            'b1rms',6.8};

        options = struct(); 

    end 


end 