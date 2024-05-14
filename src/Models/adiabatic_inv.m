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
            %  options= struct();
            % previousOptions = struct();

        
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


            % Set options as struct so when you call it a structure is created
            % for each button option 
            options= struct();
            previousOptions = struct();

        end 

        methods 

            function obj = adiabatic_inv() 
                obj.options = button2opts(obj.buttons);
                obj.previousOptions = obj.options;
                %obj = UpdateFields(obj);                
            end


            function checkfields = checkupdatedfields(obj)
                % C.R. split in two
                % checkfields = ~isequal(obj.options.TissueType, obj.previousOptions.TissueType) || ... 
                %            ~isequal(obj.options.B0, obj.previousOptions.B0) || ...
                %            ~isequal(obj.options.Pulse, obj.previousOptions.Pulse)||...
                %            obj.options.PlotAdiabatic ~= obj.previousOptions.PlotAdiabatic||... % does not equal
                %            obj.options.BlochSim1Pool ~= obj.previousOptions.BlochSim1Pool ||...
                %            obj.options.BlochSim2Pool ~= obj.previousOptions.BlochSim2Pool;
                
                % This needs to be separated otherwise it will always reset
                % params.
                
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
                    
                    % Call plotOptions function to connect changing fields 
                    %plotOptions(obj,Params);

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
                    Params.Trf = obj.Prot.PulseParameters.Mat(6);         % Trf
                    Params.nSamples = obj.Prot.PulseParameters.Mat(4);    % nSamples
                    Params.PulseOpt.beta = obj.Prot.PulseParameters.Mat(1);
                    Params.PulseOpt.A0 = obj.Prot.PulseParameters.Mat(2);
                    Params.PulseOpt.n = obj.Prot.PulseParameters.Mat(3);
                    Params.PulseOpt.Q = obj.Prot.PulseParameters.Mat(5);
                    Params.shape = obj.options.Pulse; 
                    
                    %disp(Params.A0)

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
 

% 672; 1; 5; 13.726; 512; 10.24/1000
% 1; 35; 35e-3; 0.25; 11.1e-6; 1; 0.155

             %  function setupEventListeners(obj)
             %    % Add listeners for each button press
             %    obj.buttons{4,2}.Callback = @(~,~) obj.onButtonPress('PlotAdiabatic');
             %    obj.buttons{5,2}.Callback = @(~,~) obj.onButtonPress('BlochSim1Pool');
             %    obj.buttons{6,2}.Callback = @(~,~) obj.onButtonPress('BlochSim2Pool');
             %  end
             % 
             %  function onButtonPress(obj, buttonName)
             %    % This function is called when a button is pressed
             %    %obj.UpdateFields(); % Update fields first to get the latest options
             %    obj.plotOptions(buttonName);
             %    switch buttonName
             %        case 'PlotAdiabatic'
             %            obj.plotOptions('PlotAdiabatic');
             %        case 'BlochSim1Pool'
             %            obj.plotOptions('BlochSim1Pool');
             %        case 'BlochSim2Pool'
             %            obj.plotOptions('BlochSim2Pool');
             %        otherwise
             %            error('Unknown button press');
             %    end
             %  end
             % 
             %  function plotOptions(obj,action)
             %    %plotType = obj.options;
             %    Params.Trf = obj.Prot.PulseParameters.Mat(6);         % Trf
             %    Params.nSamples = obj.Prot.PulseParameters.Mat(4);    % nSamples
             %    Params.shape = obj.options.Pulse;                     % Pulse 
             % 
             %    % Call getAdaiabatic for case to pulse
             %    [inv_pulse, omega1, A_t, Params] = getAdiabaticPulse( Params.Trf, Params.shape, Params);
             %    t = linspace(0, Params.Trf, Params.nSamples); 
             % 
             %    switch action
             % % If selecting PlotAdiabatic, call these functions and params
             %        case 'PlotAdiabatic'
             %        %Params.shape = obj.options.Pulse;                 % Pulse
             %        plotAdiabaticPulse(t, inv_pulse, A_t, omega1, Params);
             %        %obj.options.PlotAdiabatic = false;
             % 
             % % If selecting BlochSim1Pool, call these functions and params
             %        case 'BlochSim1Pool'
             %        Params.NumPools = 1;
             %        Params.M0a = obj.Prot.DefaultTissueParams.Mat(1); % M0a 
             %        Params.T2a = obj.Prot.DefaultTissueParams.Mat(3); % T2a
             %        Params.Ra = obj.Prot.DefaultTissueParams.Mat(6);  % Ra
             % 
             %        %Params.shape = obj.options.Pulse;                 % Pulse
             %        blochSimCallFunction(inv_pulse, Params)
             %        %obj.options.BlochSim1Pool = false;
             % 
             % % If selecting BlochSim2Pool, call these functions and params
             %        case 'BlochSim2Pool'
             %        Params.NumPools = 2;
             %        Params.M0a = obj.Prot.DefaultTissueParams.Mat(1); % M0a
             %        Params.R = obj.Prot.DefaultTissueParams.Mat(2);   % R
             %        Params.T2a = obj.Prot.DefaultTissueParams.Mat(3); % T2a
             %        Params.R1b = obj.Prot.DefaultTissueParams.Mat(4); % R1b
             %        Params.T2b = obj.Prot.DefaultTissueParams.Mat(5); % T2b
             %        Params.Ra = obj.Prot.DefaultTissueParams.Mat(6);  % Ra
             %        Params.M0b = obj.Prot.DefaultTissueParams.Mat(7); % M0b 
             % 
             %        %Params.shape = obj.options.Pulse;                 % Pulse 
             %        blochSimCallFunction(inv_pulse, Params)
             %        %obj.options.BlochSim2Pool = false;
             %    end 
             % end


                % obj.options.B0 = uicontrol('Style', 'popupmenu', 'String', {'3', '7', '1.5'}, ...
                % 'Tag', 'B0', 'Callback', @(src, event) obj.dropdownCallback(src, event));
                % 
                % obj.options.TissueType = uicontrol('Style', 'popupmenu', 'String', {'WM', 'GM'}, ...
                % 'Tag', 'TissueType', 'Callback', @(src, event) obj.dropdownCallback(src, event));
                % 
                % obj.options.Pulse = uicontrol('Style', 'popupmenu', 'String', {'Hs1', 'Lorentz', 'Gaussian', 'Hanning', 'Hsn', 'Sin40'}, ...
                % 'Tag', 'Pulse', 'Callback', @(src, event) obj.dropdownCallback(src, event));


            % PrevB0 = '3';
            % PrevTissueType = 'WM';
            % PrevPulse = 'Hs1';

     %  function dropdownCallback(obj,src,event) 
            %     ddValue = src.Value ;
            % 
            %     if strcmp(src.Tag, 'B0')
            %         if ~strcmp(ddValue, obj.PrevB0)
            %             obj.PrevB0 = ddValue;
            %             obj.UpdateFields();
            %         end 
            %     elseif strcmp(src.Tag, 'TissueType')
            %         if ~strcmp(ddValue, obj.PrevTissueType)
            %             obj.prevTissueType = ddValue;
            %             obj.UpdateFields();
            %         end
            %     elseif strcmp(src.Tag, 'Pulse')
            %         if ~strcmp(ddValue, obj.PrevPulse)
            %             obj.PrevPulse = ddValue;
            %             obj.UpdateFields();
            %         end 
            %     end 
            % end 

            % 
            %             function tissuetype = TissueTypedd(obj) 
            %     tissuetype = obj.options.TissueType;
            % end 
            % 
            % function b0 = B0dd(obj)
            %     b0 = obj.options.B0;
            % end 
            % 
            % function pulse = Pulsedd(obj)
            %     pulse = obj.options.Pulse;
            % end 

            % function obj = ddCallbacks(obj)
            %     set(obj.options.TissueType, 'Callback', @(src, event) obj.ddChange(src, event));
            %     set(obj.options.B0, 'Callback', @(src, event) obj.ddChange(src, event));
            %     set(obj.options.Pulse, 'Callback', @(src, event) obj.ddChange(src, event));
            % end 
            % 
            % function ddChange(obj, src, ~)
            %     switch src
            %         case obj.options.TissueType
            %             obj.options.TissueType = get(src, 'Value');
            %         case obj.options.B0 
            %             obj.options.B0 = get(src, 'Value');
            %         case obj.options.Pulse
            %             obj.options.Pulse = get(src,'Value');
            %     end 
            % 
            %     obj = UpdateFields(obj);
            % end             


            % 

            % 
            %         % Read options
            %     Params.B0 = str2double(obj.options.B0);
            %     Params.TissueType = obj.options.TissueType;
            % 
            %     % Create flags to determine which fields need updating
            %     updatePulseParams = false;
            %     updateTissueParams = false;
            % 
            %     % Update default tissue parameters if TissueType or B0 has changed
            %     if isfield(obj.options, 'TissueType') || isfield(obj.options, 'B0')
            %         updateTissueParams = true;
            %         Params = AI_defaultTissueParams(Params);
            %     end
            % 
            %     % Update pulse parameters if Pulse option has changed
            %     if isfield(obj.options, 'Pulse')
            %         updatePulseParams = true;
            %         PulseOpt = pulseparams(obj);
            %     end
            % 
            %     % Update the Prot struct based on the flags
            %     if updateTissueParams
            %         obj.Prot.DefaultTissueParams.Mat = [Params.M0a, Params.R, Params.T2a, Params.R1b, Params.T2b, Params.Ra, Params.M0b]';
            %     end
            % 
            %     if updatePulseParams
            %         obj.Prot.PulseParameters.Mat = [PulseOpt.beta, PulseOpt.A0, PulseOpt.n, PulseOpt.nSamples, PulseOpt.Q, PulseOpt.Trf]';
            %     end
            % 
            %     % Call plotOptions if necessary
            %     if any([obj.options.PlotAdiabatic, obj.options.BlochSim1Pool, obj.options.BlochSim2Pool])
            %         plotOptions(obj);
            %     end
            % 
            % end
            % 







