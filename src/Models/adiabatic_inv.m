classdef adiabatic_inv < AbstractModel 
        properties
        MRIinputs = {}; % No data needs to be downloaded
        xnames = {}; % Box names for Fitting section which I am not using  
        voxelwise = 0; % No voxel wise fitting 

        
        % Creates sections in protocol boxes.
        Prot = struct('PulseParameters', struct('Format',{{'beta(rad/s)' ; 'A0'; 'n' ;'nSamples' ;'Q' ;'Trf'}} ...
            ,'Mat',[]), ...
            'DefaultTissueParams', struct('Format',{{'M0a'; 'Raobs'; 'R'; 'T2a'; 'T1D';  'R1b'; 'T2b'; 'Ra'; 'M0b'}},'Mat', []));
     

        buttons = {'TissueType', {'WM', 'GM'},...
            'B0', {'3', '7', '1.5'}, ...
            'Pulse', {'Hs1', 'Lorentz', 'Gaussian', 'Hanning', 'Hsn', 'Sin40'},...
            'PlotAdiabatic', 'pushbutton' ...
            'BlochSim1Pool', 'pushbutton', ... 
            'BlochSim2Pool', 'pushbutton',...
            };
        options= struct();
        TissueParams = [];

        end 


        methods 

            function obj = adiabatic_inv()
            obj.options = button2opts(obj.buttons);
            
            obj = UpdateFields(obj);
            obj = plotOptions(obj);

            end

            function obj = UpdateFields(obj)
               % Call B0 and Tissue type value first 
                B0_val = str2double(obj.options.B0);
                tissuetype_val = obj.options.TissueType;
            
                % Create struct to fill tissue params 
                Params = struct; 
                Params.B0 = B0_val;
                Params.TissueType = tissuetype_val;

                Params = AI_defaultTissueParams(Params);

                % Fill default params into the object container
                 obj.Prot.DefaultTissueParams.Mat = [Params.M0a, Params.Raobs, Params.R, Params.T2a, Params.T1D, Params.R1b, Params.T2b, Params.Ra, Params.M0b]';

                 % Pulse Params
                 PulseOpt = obj.pulseparams();
                 obj.Prot.PulseParameters.Mat = [PulseOpt.beta, PulseOpt.A0, PulseOpt.n, PulseOpt.nSamples, PulseOpt.Q, PulseOpt.Trf]';

                 plotOptions(obj);
                
            end


            function PulseParams = pulseparams(obj)
            PulseOpt = struct;
            pulseType = obj.options.Pulse;
            % Creating names for each pulse to call the params associated with
            % dropdown and object containers 
                switch pulseType
                    case 'Hs1'
                        PulseOpt = obj.options;
                        PulseParams = AI_defaultHs1Params(PulseOpt);
                    case 'Lorentz'
                        PulseOpt = obj.options;
                        PulseParams = AI_defaultLorentzParams(PulseOpt);
                    case 'Gaussian'
                        PulseOpt = obj.options;
                        PulseParams = AI_defaultGaussParams(PulseOpt);
                    case 'Hanning'
                        PulseOpt = obj.options;
                        PulseParams = AI_defaultHanningParams(PulseOpt);
                    case 'Hsn'
                        PulseOpt = obj.options;
                        PulseParams = AI_defaultHsnParams(PulseOpt);
                    case 'Sin40'
                        PulseOpt = obj.options;
                        PulseParams = AI_defaultSin40Params(PulseOpt);
                    otherwise
                        error('Unknown pulse type selected');
                end
             end

          % Function to call plotting options when user presses pushbutton
             function obj = plotOptions(obj)
                Params = struct;
                Params.Trf = obj.Prot.PulseParameters.Mat(6);
                Params.nSamples = obj.Prot.PulseParameters.Mat(4);
                Params.shape = obj.options.Pulse;
                [inv_pulse, omega1, A_t, Params] = getAdiabaticPulse( Params.Trf, Params.shape, Params);
                t = linspace(0, Params.Trf, Params.nSamples); 
             % If selecting PlotAdiabatic, call these functions and params
                if obj.options.PlotAdiabatic
                    Params.shape = obj.options.Pulse;
                    plotAdiabaticPulse(t, inv_pulse, A_t, omega1, Params);
             % If selecting BlochSim1Pool, call these functions and params
                elseif obj.options.BlochSim1Pool
                    Params.NumPools = 1;
                    Params.M0a = obj.Prot.DefaultTissueParams.Mat(1);
                    Params.Ra = obj.Prot.DefaultTissueParams.Mat(8);
                    Params.T2a = obj.Prot.DefaultTissueParams.Mat(4);
                    Params.shape = obj.options.Pulse;
                    blochSimCallFunction(inv_pulse, Params)
             % If selecting BlochSim2Pool, call these functions and params
                elseif obj.options.BlochSim2Pool
                    Params.NumPools = 2;
                    Params.M0a = obj.Prot.DefaultTissueParams.Mat(1);
                    Params.M0b = obj.Prot.DefaultTissueParams.Mat(9);
                    Params.Ra = obj.Prot.DefaultTissueParams.Mat(8);
                    Params.T2a = obj.Prot.DefaultTissueParams.Mat(4);
                    Params.T2b = obj.Prot.DefaultTissueParams.Mat(7);
                    Params.R1b = obj.Prot.DefaultTissueParams.Mat(6);
                    Params.R = obj.Prot.DefaultTissueParams.Mat(3);
                    Params.shape = obj.options.Pulse;
                    blochSimCallFunction(inv_pulse, Params)
                end 
            end 

        end  
end 
 


















