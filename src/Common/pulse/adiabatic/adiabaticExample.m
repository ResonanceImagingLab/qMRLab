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

clear all

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter 
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % water pool 

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 672; % default 672 rad/s
Params.Inv.PulseOpt.n = 1; % default 1
Params.Inv.PulseOpt.mu = 5; % default 5
Params.Inv.PulseOpt.A0 = 13.726; % default 13.726 micro Tesla
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 10.24/1000; % default 10.24 ms
Params.Inv.shape = 'Hs1';


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

clear all

% Define Initial Parameters
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % White matter
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % water pool 

% Define Lorentz Parameters
Params.Inv.PulseOpt.beta = 850;  % default 850 rad/s 
Params.Inv.PulseOpt.A0 = 18; % default 18 micro Tesla
Params.Inv.PulseOpt.Q = 1e-4; % default 1e-4
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 20/1000; % default 20 ms
    % Raising Trf made omega1 look better 
Params.Inv.shape = 'Lorentz';

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

clear all

% Define Initial Parameters
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % water pool 

% Define Gaussian Parameters
Params.Inv.PulseOpt.beta = 550; % default 550 rad/s
Params.Inv.PulseOpt.A0 = 14; % default 14 micro Tesla
Params.Inv.PulseOpt.Q = 1e-4; % default 1e-4
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 10/1000; % default 10 ms
Params.Inv.shape = 'Gauss';

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
%% Hanning 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% Define Initial Parameters
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter 
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % water pool 

% Define Hanning Parameters
Params.Inv.PulseOpt.beta = 175; % default 175 rad/s
Params.Inv.PulseOpt.A0 = 14; % default 14 micro Tesla 
Params.Inv.PulseOpt.Q = 4.2e-4; % default 4.2e-4 
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 10/1000; % default 10 ms
Params.Inv.shape = 'Hanning';

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

clear all

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter 
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % water pool 

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 265; % default 265 rad/s
Params.Inv.PulseOpt.n = 8; % default 8
Params.Inv.PulseOpt.Q = 3.9e-4; % default 4e-4
Params.Inv.PulseOpt.A0 = 11; % default 10 micro Tesla
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 10/1000; % default 10.24 ms
Params.Inv.shape = 'Hsn';

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

clear all

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter 
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms
Params.NumPools = 1; % water pool 

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 200; % default 200 rad/s
Params.Inv.PulseOpt.n = 40; % default 40
Params.Inv.PulseOpt.A0 = 12; % default 12 micro Tesla
Params.Inv.PulseOpt.Q = 6.25e-7; % default 6.3e-7
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 10/1000; % default 10 ms
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

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);

%-------------------------------
%% WATER AND MACROMOLECULAR POOL %% 
%-------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hyperbolic Secant 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter 
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms (longitudinal relaxation rate water pool)
Params.R1b = 1; % 1000ms (longitudinal relaxation rate macromolecular pool)
Params.NumPools = 2; % macromolecular pool 

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 672; % default 672 rad/s
Params.Inv.PulseOpt.n = 1; % default 1 
Params.Inv.PulseOpt.mu = 5; % default 5
Params.Inv.PulseOpt.A0 = 13.726; % default 13.726 micro Tesla
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 10.24/1000; % default 10.24 ms
Params.Inv.shape = 'Hs1';


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

clear all 

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms (longitudinal relaxation rate water pool)
Params.R1b = 1; % 1000ms (longitudinal relaxation rate macromolecular pool)
Params.NumPools = 2; % macromolecular pool

% Define Lorentz Parameters
Params.Inv.PulseOpt.beta = 850;  % default 850 rad/s 
Params.Inv.PulseOpt.A0 = 18; % default 18 micro Tesla
Params.Inv.PulseOpt.Q = 1e-4; % default 1e-4
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 20/1000; % default 20 ms
    % Raising Trf made omega1 look better 
Params.Inv.shape = 'Lorentz';

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

clear all 

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms (longitudinal relaxation rate water pool)
Params.R1b = 1; % 1000ms (longitudinal relaxation rate macromolecular pool)
Params.NumPools = 2; % macromolecular pool 

% Define Gaussian Parameters
Params.Inv.PulseOpt.beta = 550; % default 550 rad/s
Params.Inv.PulseOpt.A0 = 14; % default 14 micro Tesla
Params.Inv.PulseOpt.Q = 1e-4; % default 1e-4
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 10/1000; % default 10 ms
Params.Inv.shape = 'Gauss';

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
%% Hanning 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter 
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms (longitudinal relaxation rate water pool) 
Params.R1b = 1; % 1000ms (longitudinal relaxation rate macromolecular pool)
Params.NumPools = 2; % macromolecular pool 

% Define Hanning Parameters
Params.Inv.PulseOpt.beta = 175; % default 200 rad/s
Params.Inv.PulseOpt.A0 = 14; % default 14 micro Tesla 
Params.Inv.PulseOpt.Q = 4.2e-4; % default 4.2e-4 
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 10/1000; % default 10 ms
Params.Inv.shape = 'Hanning';

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

clear all

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms (longitudinal relaxation rate water pool)
Params.R1b = 1; % 1000ms (longitudinal relaxation rate macromolecular pool) 
Params.NumPools = 2; % macromolecular pool 

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 265; % default 265 rad/s
Params.Inv.PulseOpt.n = 8; % default 8
Params.Inv.PulseOpt.Q = 3.9e-4; % default 4e-4
Params.Inv.PulseOpt.A0 = 11; % default 10 micro Tesla
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 10/1000; % default 10.24 ms
Params.Inv.shape = 'Hsn';

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

clear all

% Define Initial Parameters 
Params.B0 = 3; % Tesla
Params.TissueType = 'WM'; % white matter 
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1; % 1000ms (longitudinal relaxation rate water pool)
Params.R1b = 1; % 1000ms (longitudinal relaxation rate macromolecular pool)
Params.NumPools = 2; % macromolecular pool 

% Define Hyperbolic Secant Parameter
Params.Inv.PulseOpt.beta = 200; % default 200 rad/s
Params.Inv.PulseOpt.n = 40; % default 40
Params.Inv.PulseOpt.A0 = 12; % default 12 micro Tesla
Params.Inv.PulseOpt.Q = 6.25e-7; % default 6.3e-7
Params.Inv.nSamples = 512; % default 512
Params.Inv.Trf = 10/1000; % default 10 ms
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

% What happens in an RF pulse 
%BlochSimCallFunction(abs(inv_pulse), t, A_t, 0, Params);







