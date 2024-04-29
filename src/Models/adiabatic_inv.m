classdef adiabatic_inv < AbstractModel 
        properties
        MRIinputs = {}; % No data needs to be downloaded
        xnames = {}; % Box names for Fitting section which I am not using  
        voxelwise = 0; % No voxel wise fitting 

        %X = AI_defaultTissueParams;
        % Creates sections in protocol boxes. Fill default vals later 
        Prot = struct('PulseParams', struct('Format',{{'beta(rad/s)' ; 'A0'; 'n' ;'nSamples' ;'Q' ;'Trf'}} ...
            ,'Mat',[]), ...
            'DefaultTissueParams', struct('Format',{{'M0a'; 'Raobs'; 'R'; 'T2a'; 'T1D'; 'lineshape'; 'R1b'; 'T2b'; 'Ra'; 'M0b'; 'D'}},'Mat', []));

        % Prot = struct('DefaultTissueParams', struct('Format',{{'B0'; 'TissueType'; 'M0a'; 'Raobs'; 'R'; 'T2a'; 'T1D'; 'lineshape'; 'R1b'; 'T2b'; 'Ra'; 'M0b'; 'D'}},'Mat', []));      

        buttons = {'TissueType', {'GM', 'WM'},...
            'B0', {'3', '7', '1.5'}, ...
            'NumPools (1 or 2)', 1, ...
            'PANEL','PlotAdiabatic',1, 'Yes',false,...
            'PANEL', 'BlochSim 1 Pool', 1, 'Yes', true, ... 
            'PANEL', 'BlochSim 2 pool', 1, 'Yes', false};
        options= struct();
        TissueParams = [];
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
    % B0 = str2double(obj.options.B0);
    % TissueType = obj.options.TissueType;
    % Params = struct();
    % Params.B0 = B0;
    % Params.TissueType = TissueType;
    % Params = AI_defaultTissueParams(Params);
    % obj.Prot.DefaultTissueParams.Mat = [TissueType, num2str(B0), num2str(Params.M0a), num2str(Params.Raobs), num2str(Params.R), ...
    % num2str(Params.T2a), num2str(Params.T1D), Params.lineshape, num2str(Params.R1b), num2str(Params.T2b), num2str(Params.Ra),...
    % num2str(Params.M0b), num2str(Params.D)];


end

end 
end 



            % paramNames = {'B0'; 'TissueType'; 'M0a'; 'Raobs'; 'R'; 'T2a'; 'T1D'; 'lineshape'; 'R1b'; 'T2b'; 'Ra'; 'M0b'; 'D'};
            % paramValues = { 3 ;    {'GM'};      [];     [];    [];   [];    [];       [];      [];     [];     [];  [];   []};
            % Params = table(paramValues{:}, 'VariableNames', paramNames);
            % Params = table({'B0'}, 3, {'TissueType'}, 'GM', {'M0a'}, [], {'Raobs'}, [], {'R'}, [], {'T2a'}, [], {'T1D'}, [], {'lineshape'}, ...
            %     [], {'R1b'}, [], {'T2b'}, [], {'Ra'}, [], {'M0b'}, [], {'D'}, []);
            % Params.B0 = 3;
            % Params.TissueType = 'GM';
            % Params = AI_defaultTissueParams(Params); 
            % obj.Prot.DefaultTissueParams.Mat = Params;

       % function tissueParams = getDefaultTissueParams(obj, paramsStruct)
%     AI_defaultTissueParams(paramsStruct);
%     tissueParams = obj.Prot.DefaultTissueParams.Mat(paramsStruct);
% end 