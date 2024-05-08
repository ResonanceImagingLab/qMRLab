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
            %voxelwise = 0; % No voxel wise fitting 

        
            % Creates sections in protocol boxes: PulseParams and TissueParams
            Prot = struct('PulseParameters', struct('Format',{{'beta(rad/s)' ; 'A0'; 'n' ;'nSamples' ;'Q' ;'Trf'}} ...
                ,'Mat',[]), ...
                'DefaultTissueParams', struct('Format',{{'M0a'; 'R'; 'T2a'; 'R1b'; 'T2b'; 'Ra'; 'M0b'}},'Mat', []));

        % Prot = struct('DefaultTissueParams', struct('Format',{{'M0a'; 'R'; 'T2a'; 'R1b'; 'T2b'; 'Ra'; 'M0b'}},'Mat', []), ... 
        %     'PulseParameters', struct('Format',{{'beta(rad/s)' ; 'A0'; 'n' ;'nSamples' ;'Q' ;'Trf'}} ...
        %     ,'Mat',[]));

    
            % Creating drop box options and push buttons in Options section 
            buttons = {'TissueType', {'WM', 'GM'},...
                'B0', {'3', '7', '1.5'}, ...
                'Pulse', {'Hs1', 'Lorentz', 'Gaussian', 'Hanning', 'Hsn', 'Sin40'},...
                'PlotAdiabatic', 'pushbutton' ...
                'BlochSim1Pool', 'pushbutton', ... 
                'BlochSim2Pool', 'pushbutton',...
                };

            % Set options as struct so when you call it a structure is created
            % for each button option 
            options= struct();

        end 


        methods 

            function obj = adiabatic_inv
                obj.options = button2opts(obj.buttons);
                obj = UpdateFields(obj);
                obj = plotOptions(obj);
            end

            function obj = UpdateFields(obj)            
                % Set B0 and tissue type to options of the associated
                % dropdown 
                Params.B0 = str2double(obj.options.B0);
                Params.TissueType = obj.options.TissueType;

                % Fill in the associated Tissue params based on B0 and
                % tissue type 
                Params = AI_defaultTissueParams(Params);

                % Fill default tissue params into the object container
                obj.Prot.DefaultTissueParams.Mat = [Params.M0a, Params.R, Params.T2a, Params.R1b, Params.T2b, Params.Ra, Params.M0b]';

                % Set up Pulse Params into object container 
                PulseOpt = pulseparams(obj);
                obj.Prot.PulseParameters.Mat = [PulseOpt.beta, PulseOpt.A0, PulseOpt.n, PulseOpt.nSamples, PulseOpt.Q, PulseOpt.Trf]';

                % Call plotOptions function to connect changing fields 
                plotOptions(obj);
                
            end
                

            function PulseParams = pulseparams(obj)
                pulseType = obj.options.Pulse; % set case name to pulse option dropdown 

            % Creating names for each pulse to call the params associated with
            % dropdown and object containers 
                switch pulseType
                    case 'Hs1'
                        PulseParams = AI_defaultHs1Params(obj.options);
                    case 'Lorentz'
                        PulseParams = AI_defaultLorentzParams(obj.options);
                    case 'Gaussian'
                        PulseParams = AI_defaultGaussParams(obj.options);
                    case 'Hanning'
                        PulseParams = AI_defaultHanningParams(obj.options);
                    case 'Hsn'
                        PulseParams = AI_defaultHsnParams(obj.options);
                    case 'Sin40'
                        PulseParams = AI_defaultSin40Params(obj.options);
                    otherwise
                        error('Unknown pulse type selected');
                end
             end

          % Function to call plotting options when user presses pushbutton
          % --> Beginning set up similar to that of adiabaticExample.m
             function obj = plotOptions(obj)
                Params.Trf = obj.Prot.PulseParameters.Mat(6);         % Trf
                Params.nSamples = obj.Prot.PulseParameters.Mat(4);    % nSamples
                Params.shape = obj.options.Pulse;                     % Pulse 

                % Call getAdaiabatic for case to pulse
                [inv_pulse, omega1, A_t, Params] = getAdiabaticPulse( Params.Trf, Params.shape, Params);
                t = linspace(0, Params.Trf, Params.nSamples); 

             % If selecting PlotAdiabatic, call these functions and params
                if obj.options.PlotAdiabatic
                    Params.shape = obj.options.Pulse;                 % Pulse
                    plotAdiabaticPulse(t, inv_pulse, A_t, omega1, Params);

             % If selecting BlochSim1Pool, call these functions and params
                elseif obj.options.BlochSim1Pool
                    Params.NumPools = 1;
                    Params.M0a = obj.Prot.DefaultTissueParams.Mat(1); % M0a 
                    Params.T2a = obj.Prot.DefaultTissueParams.Mat(3); % T2a
                    Params.Ra = obj.Prot.DefaultTissueParams.Mat(6);  % Ra
                    
                    Params.shape = obj.options.Pulse;                 % Pulse
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
    
                    Params.shape = obj.options.Pulse;                 % Pulse 
                    blochSimCallFunction(inv_pulse, Params)
                end 
            end 

        end  
end 
 


















