function [rf_pulse, omega1, A_t, Params] = Sin40_pulse(Trf, Params)

%   Sin40 (n=40) Adiabatic Sin40 RF pulse function.
%   pulse = Sin40_pulse(t, Trf, PulseOpt)
%
%   B1(t) = A(t) * exp( -1i *integral(omega1(t')) dt' )
%   where A(t) is the envelope, omega1 is the frequency sweep
%
%   Phase modulation is found from taking the integral of omega1(t)
%   Frequency modulation is time derivative of phi(t)
%
%   For the case of a Sin40 pulse:
%   A(t) = A0 * (1-|sin^n(Beta*pi*t/2)|)
%   lambda = A_0^2/(Beta*Q)
%   omega1(t) = t - integral (sin^n(Beta*pi*n/2)*(1+cos^2(Beta*pi*t/2)) dt
%
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
%   PulseOpt.n: time modulation - Typical 4 for non-selective, 1 for slab
%   Reference: Matt A. Bernstein, Kevin F. Kink and Xiaohong Joe Zhou.
%              Handbook of MRI Pulse Sequences, pp. 110, Eq. 4.10, (2004)
%
%              Tannús, A. and M. Garwood (1997). "Adiabatic pulses." 
%              NMR in Biomedicine 10(8): 423-434.
%                   --> A(t), omega1 
%
%              Kupce, E. and Freeman, R (1995). "Optimized Adiabatic Pulses
%              for Wideband Spin Inversion." Journal of Magnetic Resonance
%              Imaging, Series A 118(2): 299-303.
%                   --> lambda
%
%              Tannus, A. Garwood, M. (1996). "Improved Performance of 
%              Frequency Swept Pulses Using Offset-Independent
%              Adiabaticity" Journal of Magnetic Resonance, 120(1),
%              133-137.
%                   --> Fig 1a and 1b. Show how width of amplitude and
%                   frequency vary with each pulse
%
%              De Graaf, R. Nicolay, K. (1997). "Adiabatic rf Pulses: 
%              Applications to In Vivo NMR" Concepts Magn. Reson. 9:247-268.
%                   --> Fig 12 b. Width of inversion profile is wider than
%                   for that of hyperbolic secant 
%
%
% To be used with qMRlab
% Written by Christopher Rowley 2023 and Amie Demmans 2024


% Function to fill default values;
Params.PulseOpt = defaultSin40Params(Params.PulseOpt);

nSamples = Params.PulseOpt.nSamples;  
t = linspace(0, Trf, nSamples);

%tau = ((2*t/Trf)-1);
tau = t-Trf/2;

% Amplitude
At1 = sin(Params.PulseOpt.beta*pi.*tau./2).^Params.PulseOpt.n;
At2 = 1-abs(At1);
A_t =  Params.PulseOpt.A0* At2;
A_t((t < 0 | t>Trf)) = 0;
% disp( ['Average B1 of the pulse is:', num2str(mean(A_t))]) 

% Scaling Factor 
lambda = (Params.PulseOpt.A0)^2 ./ (Params.PulseOpt.beta.*Params.PulseOpt.Q);

% Frequency modulation function 
% Carrier frequency modulation function w(t):
omegaterm1 = At1;
omegaterm2 = 1 + cos(Params.PulseOpt.beta*pi.*tau./2).^2;
omegaterm3 = omegaterm1.*omegaterm2;
omegaterm4 = -lambda*(tau - cumtrapz(tau,omegaterm3));
omega1 = omegaterm4 - omegaterm4(round(nSamples/2));
omega = lambda*(tau-cumtrapz(tau,omegaterm3));

% Phase modulation function phi(t):
phi = cumtrapz(tau, omega);


% Put together complex RF pulse waveform:
rf_pulse = A_t .* exp(1i .* phi);










