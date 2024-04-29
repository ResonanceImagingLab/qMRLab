classdef adiabatic_inv < AbstractModel 
        properties
        MRIinputs = {}; % No data needs to be downloaded
        xnames = {}; % Box names for Fitting section which I am not using  
        voxelwise = 0; % No voxel wise fitting 

        %X = AI_defaultTissueParams;
        % Creates sections in protocol boxes. Fill default vals later 
        Prot = struct('PulseParameters', struct('Format',{{'beta(rad/s)' ; 'A0'; 'n' ;'nSamples' ;'Q' ;'Trf'}} ...
            ,'Mat',[]), ...
            'DefaultTissueParams', struct('Format',{{'M0a'; 'Raobs'; 'R'; 'T2a'; 'T1D';  'R1b'; 'T2b'; 'Ra'; 'M0b'}},'Mat', []));
     

        buttons = {'TissueType', {'WM', 'GM'},...
            'B0', {'3', '7', '1.5'}, ...
            'Pulse', {'Hs1', 'Lorentz', 'Gaussian', 'Hanning', 'Hsn', 'Sin40'},...
            'NumPools (1 or 2)', 1, ...
            'PANEL','PlotAdiabatic',1, 'Yes',false,...
            'PANEL', 'BlochSim 1 Pool', 1, 'Yes', true, ... 
            'PANEL', 'BlochSim 2 pool', 1, 'Yes', false};
        options= struct();
        TissueParams = [];
        PulseParams = [];
        end 


methods (Hidden=true)
% Hidden methods goes here.
end


methods 

function obj = adiabatic_inv()
            obj.options = button2opts(obj.buttons);
            obj = UpdateFields(obj);

end

function obj = UpdateFields(obj)
    obj.TissueParams = AI_defaultTissueParams(str2double(obj.options.B0),obj.options.TissueType);
    obj.Prot.DefaultTissueParams.Mat = obj.TissueParams;
    obj.PulseParams = obj.pulseparams();
    obj.Prot.PulseParameters.Mat = obj.PulseParams;

end

function PulseParams = pulseparams(obj)
    pulseType = obj.options.Pulse;
    switch pulseType
        case 'Hs1'
            PulseParams = AI_defaultHs1Params();
        case 'Lorentz'
            PulseParams = AI_defaultLorentzParams();
        case 'Gaussian'
            PulseParams = AI_defaultGaussParams();
        case 'Hanning'
            PulseParams = AI_defaultHanningParams();
        case 'Hsn'
            PulseParams = AI_defaultHsnParams();
        case 'Sin40'
            PulseParams = AI_defaultSin40Params();
        otherwise
            error('Unknown pulse type selected');
    end
end

end 
end 


















