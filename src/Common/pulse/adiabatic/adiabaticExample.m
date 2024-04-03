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
%PlotAdiabaticPulse(t, inv_pulse);

% Plot Bloch Sim Results based on NumPools
BlochSimCallFunction(inv_pulse, t, A_t, omega1,Params)

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lorentz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 1000/80; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1;

% Define Lorentz Parameters
Params.Inv.PulseOpt.beta = 600;  
Params.Inv.PulseOpt.A0 = 18;
Params.Inv.PulseOpt.Q = 1e-4;
    % 1e-4 will give nice boat shape but HUGE omega1
    % 1e-3 gives omega1 value better but Mz is crazy 
Params.Inv.nSamples = 512;
Params.Inv.Trf = 20/1000;
    % Raising Trf made omega1 look better 
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

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Initial Parameters
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 1000/80; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 2;

% Define Gaussian Parameters
Params.Inv.PulseOpt.beta = 550; 
    % B = 600 --> happy
Params.Inv.PulseOpt.A0 = 13;
    % A0 = 14 --> happy 
Params.Inv.PulseOpt.Q = 1e-4;
    % 1e-4 gives the bowl look but very high omega (3200 Hz)
    % 2e-4 gives more of a v shaped bowl 
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
    % omega1 is unitless based on equation 

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hanning 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Initial Parameters
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 1000/80; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1;

% Define Hanning Parameters
Params.Inv.PulseOpt.beta = 175;
    % Can vary from 150-200 depending on what you're looking for 
Params.Inv.PulseOpt.A0 = 14;
Params.Inv.PulseOpt.Q = 4.2e-4;
    % Any deviations from this and Mz got mad 
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000;
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

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hsn (Params set for n = 8) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 1000/80; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % or 2

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 250; 
% amplitude looks best at beta = 250 
Params.Inv.PulseOpt.n = 8; 
Params.Inv.PulseOpt.Q = 4e-4;
Params.Inv.PulseOpt.A0 = 10;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000; 
Params.Inv.shape = 'Hsn8';

% Need to figure out scaling on y-axis of Mz plot as well as issues with
% spiking on the bottom --> Is this from scaling issues? 

% Apply Inversion pulse by calling GetAdiabatic
Params.dispFigure = 0;
[inv_pulse, omega1, A_t, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% To check your pulse: Plot  
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
%PlotAdiabaticPulse(t, inv_pulse);

% Plot Bloch Sim Results based on NumPools 
BlochSimCallFunction(inv_pulse, t, A_t, omega1, Params);

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sin40 (n = 40) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % or 2

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 200; 
% DO NOT stray from 200 it gets angry!!!!
Params.Inv.PulseOpt.n = 40; 
Params.Inv.PulseOpt.A0 = 14;
Params.Inv.PulseOpt.Q = 8e-7;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000;
Params.Inv.shape = 'Sin40';


% Apply Inversion pulse by calling GetAdiabatic
Params.dispFigure = 0;
[inv_pulse, omega1, A_t, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% To check your pulse: Plot  
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
%PlotAdiabaticPulse(t, inv_pulse);

% Plot Bloch Sim Results based on NumPools 
BlochSimCallFunction(inv_pulse, t, A_t, omega1, Params);
