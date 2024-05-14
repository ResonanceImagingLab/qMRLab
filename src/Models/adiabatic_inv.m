classdef adiabatic_inv < AbstractModel 
% adiabatic_inv: Adiabatic inversion pulses 
%
% Assumptions: 
%
% Inputs: 
%
% Outputs:
%  PlotAdiabaticPulse          View the amplitude modulation, phase
%                              modulation and frequency modulation of each pulse 
%
%  BlochSim1Pool               Assess the inversion characteristics of the
%                              selected pulse using bloch simulations for a 
%                              water pool case 
%
%  BlochSim2Pool               Asses the inversion characteristics of the
%                              selected pulse using bloch simulations for 
%                              the water and bound pool case
%
% Options: 
%  TissueType                  Select the desired tissue type 
%                               - White matter (WM)
%                               - Grey matter (GM) 
%
%  B0                          Select the magnet size
%                               - 3 T
%                               - 7 T 
%                               - 1.5 T 
%  
% Pulse                        Select the pulse type you wish to view 
%                               - Hs1
%                               - Lorentz
%                               - Gaussian 
%                               - Hanning 
%                               - Hsn
%                               - Sin40
%
% Plotting Option              Select the plotting options (refer to outputs for descriptions) 
%                               - PlotAdiabatic 
%                               - BlochSim1Pool 
%                               - BlochSim2Pool 
%
% Authors: Amie Demmans, 2024 
%
% References:
%   Please refer to getAdiabaticPulse.m for all references used to develop
%     this module
%   In addition to the citing package: 
%     Karakuzu A., Boudreau M., Duval T.,Boshkovski T., Leppert I.R., Cabana J.F., 
%     Gagnon I., Beliveau P., Pike G.B., Cohen-Adad J., Stikov N. (2020), qMRLab: 
%     Quantitative MRI analysis, under one umbrella doi: 10.21105/joss.02343



properties
    MRIinputs = {}; % No data needs to be downloaded
    xnames = {}; % Box names for Fitting section which I am not using, still needs to be defined though


    % Creates sections in protocol boxes: PulseParams and TissueParams
    Prot = struct('PulseParameters', struct('Format',{{'beta(rad/s)' ; 'A0'; 'n' ;'nSamples' ;'Q' ;'Trf'}} ...
        ,'Mat',[672; 13.726; 1; 512; 5; 10.24/1000]), ...
        'DefaultTissueParams', struct('Format',{{'M0a'; 'R'; 'T2a'; 'R1b'; 'T2b'; 'Ra'; 'M0b'}}, ...
        'Mat', [1; 35; 35e-3; 0.25; 11.1e-6; 1; 0.155]));


    % Creating drop box options and push buttons in Options section
    buttons = {'TissueType', {'WM', 'GM'},...
        'B0', {'3', '7', '1.5'}, ...
        'Pulse', {'Hs1', 'Lorentz', 'Gaussian', 'Hanning', 'Hsn', 'Sin40'},...
        'PlotAdiabatic', 'pushbutton', ...
        'BlochSim1Pool', 'pushbutton', ...
        'BlochSim2Pool', 'pushbutton',...
        };


    % Set options and previousOptions as struct so when you call it
    % a structure is created for each button option
    options= struct();
    previousOptions = struct();

end

methods

    function obj = adiabatic_inv()
        obj.options = button2opts(obj.buttons);
        obj.previousOptions = obj.options;
    end


    function checkfields = checkupdatedfields(obj)
        if (~isequal(obj.options.TissueType, obj.previousOptions.TissueType) || ...
                ~isequal(obj.options.B0, obj.previousOptions.B0) || ...
                ~isequal(obj.options.Pulse, obj.previousOptions.Pulse))
            checkfields = 1; % reset params to defaults
        elseif (obj.options.PlotAdiabatic ~= obj.previousOptions.PlotAdiabatic||... % does not equal
                obj.options.BlochSim1Pool ~= obj.previousOptions.BlochSim1Pool ||...
                obj.options.BlochSim2Pool ~= obj.previousOptions.BlochSim2Pool)
            checkfields = 2; % run sims -> this needs to be moved to a new function C.R.
        else
            checkfields = 0;
        end

        % Debug C.R.
        str = ['checkfields = ', num2str(checkfields)];
        disp(str);

    end


    function obj = UpdateFields(obj)

        % Debug C.R.
        disp('Update Fields called with checkFields');
        disp(obj.checkupdatedfields)

        if obj.checkupdatedfields == 1 % C.R. add checkfields

            %obj.Storedparams = Params;
            obj.previousOptions = obj.options;

            %Set B0 and tissue type to options of the associated
            %dropdown
            Params.B0 = str2double(obj.options.B0);
            Params.TissueType = obj.options.TissueType;

            % Fill in the associated Tissue params based on B0 and
            % tissue type
            Params = AI_defaultTissueParams(Params);

            % Fill default tissue params into the object container
            obj.Prot.DefaultTissueParams.Mat = [Params.M0a, Params.R, Params.T2a, Params.R1b, Params.T2b, Params.Ra, Params.M0b]';

            % Set up Pulse Params into object container
            PulseOpt = pulseparams(obj);
            obj.Prot.PulseParameters.Mat = [PulseOpt.beta, PulseOpt.A0, PulseOpt.n, PulseOpt.nSamples, PulseOpt.Q, PulseOpt.Trf]' ;

        elseif obj.checkupdatedfields == 2
            plotOptions(obj);
        end

    end


    function obj = pulseparams(obj)
        pulseType = obj.options.Pulse; % set case name to pulse option dropdown

        % Creating names for each pulse to call the params associated with
        % dropdown and object containers
        switch pulseType
            case 'Hs1'
                obj = AI_defaultHs1Params(obj.options);
            case 'Lorentz'
                obj = AI_defaultLorentzParams(obj.options);
            case 'Gaussian'
                obj = AI_defaultGaussParams(obj.options);
            case 'Hanning'
                obj = AI_defaultHanningParams(obj.options);
            case 'Hsn'
                obj = AI_defaultHsnParams(obj.options);
            case 'Sin40'
                obj = AI_defaultSin40Params(obj.options);
            otherwise
                error('Unknown pulse type selected');
        end
    end


    %Function to call plotting options when user presses pushbutton
    %--> Beginning set up similar to that of adiabaticExample.m
    function obj = plotOptions(obj)
        % Call Trf, nSamples and shape for getAdiabatic
        Params.Trf = obj.Prot.PulseParameters.Mat(6);         % Trf
        Params.nSamples = obj.Prot.PulseParameters.Mat(4);    % nSamples
        Params.shape = obj.options.Pulse;                     % pulseshape
        % Call other pulse properties to allow editing 
        Params.PulseOpt.beta = obj.Prot.PulseParameters.Mat(1); % beta
        Params.PulseOpt.A0 = obj.Prot.PulseParameters.Mat(2);   % A0
        Params.PulseOpt.n = obj.Prot.PulseParameters.Mat(3);    % n
        Params.PulseOpt.Q = obj.Prot.PulseParameters.Mat(5);    % Q
        

        % Call getAdaiabatic for case to pulse
        [inv_pulse, omega1, A_t, Params] = getAdiabaticPulse( Params.Trf, Params.shape, Params);
        t = linspace(0, Params.Trf, Params.nSamples);


        % If selecting PlotAdiabatic, call these functions and params
        if obj.options.PlotAdiabatic
            plotAdiabaticPulse(t, inv_pulse, A_t, omega1, Params);

        % If selecting BlochSim1Pool, call these functions and params
        elseif obj.options.BlochSim1Pool
            Params.NumPools = 1;
            Params.M0a = obj.Prot.DefaultTissueParams.Mat(1); % M0a
            Params.T2a = obj.Prot.DefaultTissueParams.Mat(3); % T2a
            Params.Ra = obj.Prot.DefaultTissueParams.Mat(6);  % Ra

            blochSimCallFunction(inv_pulse, Params)

        % If selecting BlochSim2Pool, call these functions and params
        elseif obj.options.BlochSim2Pool
            Params.NumPools = 2;
            Params.M0a = obj.Prot.DefaultTissueParams.Mat(1); % M0a
            Params.R = obj.Prot.DefaultTissueParams.Mat(2);   % R
            Params.T2a = obj.Prot.DefaultTissueParams.Mat(3); % T2a
            Params.R1b = obj.Prot.DefaultTissueParams.Mat(4); % R1b
            Params.T2b = obj.Prot.DefaultTissueParams.Mat(5); % T2b
            Params.Ra = obj.Prot.DefaultTissueParams.Mat(6);  % Ra
            Params.M0b = obj.Prot.DefaultTissueParams.Mat(7); % M0b

            blochSimCallFunction(inv_pulse, Params)

        end

    end

end
end
 















