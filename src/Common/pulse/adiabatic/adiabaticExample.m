%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hyperbolic Secant 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % or 2

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 672;  
Params.Inv.PulseOpt.n = 1; 
Params.Inv.PulseOpt.mu = 5;
Params.Inv.PulseOpt.A0 = 13.726;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10.24/1000;
Params.Inv.shape = 'hsn';


% Apply Inversion pulse by calling GetAdiabatic
Params.dispFigure = 0;
[inv_pulse, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% To check your pulse: Plot  
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
PlotPulse(t, inv_pulse);

% Plot Bloch Sim Results based on NumPools
[rf_pulse, omega1, A_t, ~] = hyperbolicSecant_pulse(Params.Inv.Trf, Params.Inv);

BlochSimCallFunction(inv_pulse, t, A_t, omega1,Params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lorentz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.NumPools = 1;

% Define Lorentz Parameters
Params.Inv.PulseOpt.beta = 99;  
Params.Inv.PulseOpt.A0 = 13.726;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 15/1000;
Params.Inv.shape = 'lorentz';

% Apply inversion pulse by calling GetAdiabatic 
Params.dispFigure = 0;
[inv_pulse, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% To check your pulse: Plot  
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
PlotPulse(t, inv_pulse);

M_start = [0,0,Params.M0a]'; % starting mangetization 
delta_start = -2000;
delta_end = 2000; 
delta = linspace(delta_start,delta_end,length(t));

M_return = blochSimAdiabaticPulse_1pool(inv_pulse, Params.Inv.Trf, delta, Params, M_start);
disp('Final magnetization after inversion pulse')
disp(M_return)

% Bloch sim and return magnetization
%
% Params.PulseOpt = defaultLorentzParams(Params); % Define PulseOpt Parameters
% Params.Inv.PulseOpt = Params.PulseOpt; % Set Inverse PulseOpt Params to same as PulseOpt 
% [rf_pulse, omega1, A_t, ~] = Lorentz_pulse(Params.Inv.Trf, Params.Inv); % Call Hyperbolicsecant 
% 
% 
% BlochSimCallFunction(rf_pulse, t, A_t, omega1, Params);