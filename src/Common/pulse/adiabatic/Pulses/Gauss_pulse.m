function [rf_pulse, omega1, A_t, Params] = Gauss_pulse( Trf, Params)

%   GaussC_pulse Adiabatic Gaussian RF pulse function.
%   pulse = GaussC_pulse(Trf, PulseOpt)
%
%   B1(t) = A(t) * exp( -1i *integral(omega1(t')) dt' )
%   where A(t) is the envelope, omega1 is the frequency sweep
%
%   Phase modulation is found from taking the integral of omega1(t)
%   Frequency modulation is time derivative of phi(t)
%
%   For the case of a Gauus^c pulse:
%   A(t) = A_0 * exp(-Beta * t^2)
%   lambda = A0^2/(Beta*Q)
%   omega1(t) = lamdba*erf(Beta*t)
%   A0 is the peak amplitude in microTesla
%   Beta is a frequency modulation parameter in rad/s
%
%   The pulse is defined to be 0 outside the pulse window (before 
%   t = 0 or after t=Trf). (HSn, n = 1-8+) 
%
%   --args--
%   t: Function handle variable, represents the time.
%   Trf: Duration of the RF pulse in seconds.
%
%   --optional args--
%   PulseOpt: Struct. Contains optional parameters for pulse shapes.
%   PulseOpt.Beta: frequency modulation parameter
% 
%   Reference: Matt A. Bernstein, Kevin F. Kink and Xiaohong Joe Zhou.
%              Handbook of MRI Pulse Sequences, pp. 110, Eq. 4.10, (2004)
%
%              Tannús, A. and M. Garwood (1997). "Adiabatic pulses." 
%              NMR in Biomedicine 10(8): 423-434.
%                  --> Table 1 contains all modulation functions 
%                  --> A(t), omega1
%                  --> Fig 5, Gaussian OIA pulse image 
%                  --> Trf = 10 ms 
%
%              Kupce, E. and Freeman, R (1995). "Optimized Adiabatic Pulses
%              for Wideband Spin Inversion." Journal of Magnetic Resonance
%              Imaging, Series A 118(2): 299-303.
%                  --> lambda equation, Eq. 10 (added to omega1 for scaling)
%
%              Tannus, A. Garwood, M. (1996). "Improved Performance of 
%              Frequency Swept Pulses Using Offset-Independent
%              Adiabaticity" Journal of Magnetic Resonance, 120(1),
%              133-137.
%                   --> Fig 1a and 1b. Show how width of amplitude and
%                   frequency vary with each pulse
%
% To be used with qMRlab
% Written by Christopher Rowley 2023 & Amie Demmans 2024


% Function to fill default values;
Params.PulseOpt = defaultGaussParams(Params.PulseOpt);

nSamples = Params.PulseOpt.nSamples;  
t = linspace(0, Trf, nSamples);
tau = t-Trf/2;

% Amplitude --> From Ref 2
A_t = Params.PulseOpt.A0 .* exp(-1*(Params.PulseOpt.beta.^2 .* tau.^2)./2);
A_t((t < 0 | t>Trf)) = 0;
% disp( ['Average B1 of the pulse is:', num2str(mean(A_t))]) 

% Scaling Factor 
lambda = (Params.PulseOpt.A0)^2 ./ (Params.PulseOpt.beta.*Params.PulseOpt.Q);

% Carrier frequency modulation function w(t):  
omega1 = -lambda*erf(Params.PulseOpt.beta.*tau)./erf(Params.PulseOpt.beta);

% Phase modulation function phi(t):
phi1num = lambda.*tau.*erf(Params.PulseOpt.beta.*tau);
phi1denom = erf(Params.PulseOpt.beta);
phi1 = phi1num/phi1denom;
phi2num = lambda*exp(-Params.PulseOpt.beta.^2 .* tau.^2);
phi2denom = sqrt(pi)*Params.PulseOpt.beta*erf(Params.PulseOpt.beta); % If I divide B by 2pi for rad/s to Hz it gets angry 
phi2 = phi2num/phi2denom;
phi = phi1+phi2;

% Put together complex RF pulse waveform:
rf_pulse = A_t .* exp(1i .* phi);




