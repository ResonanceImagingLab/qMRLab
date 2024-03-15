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
[inv_pulse, omega1, A_t, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% To check your pulse: Plot  
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
PlotAdiabaticPulse(t, inv_pulse);

% Plot Bloch Sim Results based on NumPools

BlochSimCallFunction(inv_pulse, t, A_t, omega1,Params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lorentz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1;

% Define Lorentz Parameters
Params.Inv.PulseOpt.beta = 100;  
Params.Inv.PulseOpt.A0 = 11;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 15/1000;
Params.Inv.shape = 'lorentz';

% Apply inversion pulse by calling GetAdiabatic 
Params.dispFigure = 0;
[inv_pulse, omega1, A_t, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% To check your pulse: Plot  
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
%PlotAdiabaticPulse(t, inv_pulse);

% Plot Bloch Sim Results based on NumPools 
BlochSimCallFunction(inv_pulse, t, A_t, omega1, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Initial Parameters
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1;

% Define Gaussian Parameters
Params.Inv.PulseOpt.beta = 2;  
Params.Inv.PulseOpt.A0 = 15;
Params.Inv.PulseOpt.Q = 4;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000;
Params.Inv.shape = 'gauss';

% Apply inversion pulse by calling GetAdiabatic 
Params.dispFigure = 0;
[inv_pulse, omega1, A_t, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% To check your pulse: Plot  
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
%PlotAdiabaticPulse(t, inv_pulse);

% Plot Bloch Sim Results based on NumPools 
BlochSimCallFunction(inv_pulse, t, A_t, omega1, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hanning 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Initial Parameters
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1;

% Define Gaussian Parameters 
Params.Inv.PulseOpt.A0 = 6.2;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 2/1000;
Params.Inv.shape = 'hanning';

% Apply inversion pulse by calling GetAdiabatic 
Params.dispFigure = 0;
[inv_pulse, omega1, A_t, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% To check your pulse: Plot  
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
%PlotAdiabaticPulse(t, inv_pulse);

% Plot Bloch Sim Results based on NumPools 
BlochSimCallFunction(inv_pulse, t, A_t, omega1, Params);


