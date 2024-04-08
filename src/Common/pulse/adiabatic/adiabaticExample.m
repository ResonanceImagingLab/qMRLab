%% ADIABATIC INVERSION PULSES %%
% ------------------------------
% This code is set to call functions that will plot the different adiabatic
% inversion pulses.
%
% 1. Each section is set for the six different adiabatic inversion pulses 
%    - Hyperbolic Secant (HS1) 
%    - Lorentzian 
%    - Gaussian 
%    - Hanning 
%    - HSn (n = 2-8) *default params set for n = 8*
%    - Sin40 (n = 40) 
%
% 2. The pulses are listed twice either for a 1 pool or 2 pool case 
%    - First set of six are set to Params.NumPools = 1 (water only)
%    - Second set of six are set to Params.NumPools = 2 (water & bound)
%
% 3. There are three different plotting options for each pulse:
%    - Plot Adiabatic pulse to view amplitude and frequency modulation 
%    - Apply Bloch sim equations to visualize end of inversion pulse 
%      with specified B1 and B0 inhomogeneities  
%    - Plot just the RF pulse to visualize when frequency modulation is absent 
%
% 4. The exact amplitude, frequency, and phase modulation equation for each
%    pulse with their references are listed in the matlab functions
%
% 5. Pulse parameters are automatically set with the default params, but
%    you can change them as you please to see how that may affect the pulse
%    - Descriptions of each parameter are listed in the defaultparam
%    matlab functions for each pulse 
%    - Units for each parameter are listed in this code where applicable
%
% 6. For more information on adiabatic inversion pulses, check out the
%    README.md or all references are listed in GetAdiabaticPulse.m
%
% Written by Amie Demmans 2024
%
%--------------------
%% WATER POOL ONLY %% 
%--------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hyperbolic Secant 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % or 2

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 672; % rad/s
Params.Inv.PulseOpt.n = 1; 
Params.Inv.PulseOpt.mu = 5;
Params.Inv.PulseOpt.A0 = 13.726; % micro Tesla
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10.24/1000; % ms
Params.Inv.shape = 'hs1';


% Apply Inversion pulse by calling GetAdiabatic
Params.dispFigure = 0;
[inv_pulse, omega1, A_t, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% To check your pulse: Plot  
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
%PlotAdiabaticPulse(t, inv_pulse);

% Plot Bloch Sim Results based on NumPools
BlochSimCallFunction(inv_pulse, t, A_t, omega1,Params)
title('HS1')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lorentz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 1000/80; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1;

% Define Lorentz Parameters
Params.Inv.PulseOpt.beta = 600;  % rad/s 
Params.Inv.PulseOpt.A0 = 18; % micro Tesla
Params.Inv.PulseOpt.Q = 1e-4;
    % 1e-4 will give nice boat shape but HUGE omega1
    % 1e-3 gives omega1 value better but Mz is crazy 
Params.Inv.nSamples = 512;
Params.Inv.Trf = 20/1000; % ms
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
title('Lorentz')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 1000/80; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1;

% Define Gaussian Parameters
Params.Inv.PulseOpt.beta = 550; % rad/s
    % B = 600 --> happy
Params.Inv.PulseOpt.A0 = 13; % micro Tesla
    % A0 = 14 --> happy 
Params.Inv.PulseOpt.Q = 1e-4;
    % 1e-4 gives the bowl look but very high omega (3200 Hz)
    % 2e-4 gives more of a v shaped bowl 
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000; % ms
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
title('Gaussian')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hanning 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 1000/80; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1;

% Define Hanning Parameters
Params.Inv.PulseOpt.beta = 175; % rad/s
    % Can vary from 150-200 depending on what you're looking for 
Params.Inv.PulseOpt.A0 = 14; % micro Tesla 
Params.Inv.PulseOpt.Q = 4.2e-4;
    % Any deviations from this and Mz got mad 
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000; % ms
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
title('Hanning')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hsn (Params set for n = 8) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 1000/80; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % or 2

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 250; % rad/s
% amplitude looks best at beta = 250 
Params.Inv.PulseOpt.n = 8; 
Params.Inv.PulseOpt.Q = 4e-4;
Params.Inv.PulseOpt.A0 = 10; % micro Tesla
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000; % ms
Params.Inv.shape = 'Hsn';

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
title('HSn')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sin40 (n = 40) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % or 2

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 200; % rad/s
% DO NOT stray from 200 it gets angry!!!!
Params.Inv.PulseOpt.n = 40; 
Params.Inv.PulseOpt.A0 = 13; % micro Tesla
Params.Inv.PulseOpt.Q = 6.5e-7;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000; % ms
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
title('Sin40')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%--------------------------
%% WATER AND BOUND POOL %% 
%--------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hyperbolic Secant 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.R2b = 0.012; % 0.012 ms
Params.Ra = 1; % 1000ms
Params.R1b = 1; % 1000ms
Params.kr = 200*1000; % 200 /ms
Params.kf = 0.2*1000; % 2 /ms
Params.NumPools = 2; 

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 672; % rad/s
Params.Inv.PulseOpt.n = 1; 
Params.Inv.PulseOpt.mu = 5;
Params.Inv.PulseOpt.A0 = 13.726; % micro Tesla
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10.24/1000; % ms
Params.Inv.shape = 'hs1';


% Apply Inversion pulse by calling GetAdiabatic
Params.dispFigure = 0;
[inv_pulse, omega1, A_t, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% To check your pulse: Plot  
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
%PlotAdiabaticPulse(t, inv_pulse);

% Plot Bloch Sim Results based on NumPools
BlochSimCallFunction(inv_pulse, t, A_t, omega1,Params)
title('HS1')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lorentz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.R2b = 12/1e6; % 12 us
Params.Ra = 1; % 1000ms
Params.R1b = 1; % 1000ms
Params.kr = 200*1000; % 200 /ms
Params.kf = 0.2*1000; % 2 /ms
Params.NumPools = 2; 

% Define Lorentz Parameters
Params.Inv.PulseOpt.beta = 600;  % rad/s 
Params.Inv.PulseOpt.A0 = 18; % micro Tesla
Params.Inv.PulseOpt.Q = 1e-4;
    % 1e-4 will give nice boat shape but HUGE omega1
    % 1e-3 gives omega1 value better but Mz is crazy 
Params.Inv.nSamples = 512;
Params.Inv.Trf = 20/1000; % ms
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
title('Lorentz')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.R2b = 12/1e6; % 12 us
Params.Ra = 1; % 1000ms
Params.R1b = 1; % 1000ms
Params.kr = 200*1000; % 200 /ms
Params.kf = 0.2*1000; % 2 /ms
Params.NumPools = 2; 

% Define Gaussian Parameters
Params.Inv.PulseOpt.beta = 550; % rad/s
    % B = 600 --> happy
Params.Inv.PulseOpt.A0 = 13; % micro Tesla
    % A0 = 14 --> happy 
Params.Inv.PulseOpt.Q = 1e-4;
    % 1e-4 gives the bowl look but very high omega (3200 Hz)
    % 2e-4 gives more of a v shaped bowl 
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000; % ms
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
title('Gaussian')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hanning 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.R2b = 12/1e6; % 12 us
Params.Ra = 1; % 1000ms
Params.R1b = 1; % 1000ms
Params.kr = 200*1000; % 200 /ms
Params.kf = 0.2*1000; % 2 /ms
Params.NumPools = 2; 

% Define Hanning Parameters
Params.Inv.PulseOpt.beta = 175; % rad/s
    % Can vary from 150-200 depending on what you're looking for 
Params.Inv.PulseOpt.A0 = 14; % micro Tesla 
Params.Inv.PulseOpt.Q = 4.2e-4;
    % Any deviations from this and Mz got mad 
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000; % ms
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
title('Hanning')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hsn (Params set for n = 8) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.R2b = 12/1e6; % 12 us
Params.Ra = 1; % 1000ms
Params.R1b = 1; % 1000ms
Params.kr = 200*1000; % 200 /ms
Params.kf = 0.2*1000; % 2 /ms
Params.NumPools = 2; 

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 250; % rad/s
% amplitude looks best at beta = 250 
Params.Inv.PulseOpt.n = 8; 
Params.Inv.PulseOpt.Q = 4e-4;
Params.Inv.PulseOpt.A0 = 10; % micro Tesla
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000; % ms
Params.Inv.shape = 'Hsn';

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
title('HSn')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sin40 (n = 40) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.R2a = 80/1000; % 80 ms
Params.R2b = 12/1e6; % 12 us
Params.Ra = 1; % 1000ms
Params.R1b = 1; % 1000ms
Params.kr = 200*1000; % 200 /ms
Params.kf = 0.2*1000; % 2 /ms
Params.NumPools = 2; 

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 200; % rad/s
% DO NOT stray from 200 it gets angry!!!!
Params.Inv.PulseOpt.n = 40; 
Params.Inv.PulseOpt.A0 = 13; % micro Tesla
Params.Inv.PulseOpt.Q = 6.5e-7;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10/1000; % ms
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
title('Sin40')

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);
