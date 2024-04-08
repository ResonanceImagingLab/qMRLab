function [rf_pulse, omega1, A_t, Params] = Hsn_pulse(Trf, Params)

%   hyperbolicSecant_pulse Adiabatic hyperbolic secant RF pulse function.
%   pulse = hyperbolicSecant_pulse(t, Trf, PulseOpt)
%
%   B1(t) = A(t) * exp( -1i *integral(omega1(t')) dt' )
%   where A(t) is the envelope, omega1 is the frequency sweep
%
%   Phase modulation is found from taking the integral of omega1(t)
%   Frequency modulation is time derivative of phi(t)
%
%   For the case of a hyperbolic secant pulse:
%   A(t) = A0 * sech((Beta*t)^n)
%   lambda = (A_0^2/(Beta*Q))^2
%   omega1(t) = integral(sech^2(Beta*t^n)
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
%                   --> Table 1 contains all modulation functions 
%                   --> A(t), omega1 
%
%              Kupce, E. and Freeman, R (1995). "Optimized Adiabatic Pulses
%              for Wideband Spin Inversion." Journal of Magnetic Resonance
%              Imaging, Series A 118(2): 299-303.
%                   --> lambda equation added to omega 1 for scaling, Eq.10
%
%              Tesiram, Y. "Implementation Equations for HSn RF Pulses."
%              Journal of Magentic Resonance, 204, 333-339
%                   --> Placing beta ^n 
%
%              Tannus, A. Garwood, M. (1996). "Improved Performance of 
%              Frequency Swept Pulses Using Offset-Independent
%              Adiabaticity" Journal of Magnetic Resonance, 120(1),
%              133-137.
%                   --> Fig 1a and 1b. Show how width of amplitude and
%                   frequency vary with each pulse
%
%
% To be used with qMRlab
% Written by Christopher Rowley 2023 and Amie Demmans 2024


% Function to fill default values;
Params.PulseOpt = defaultHsnParams(Params.PulseOpt);

nSamples = Params.PulseOpt.nSamples;  
t = linspace(0, Trf, nSamples);
tau = t-Trf/2;

% Amplitude
A_t =  Params.PulseOpt.A0* sech((Params.PulseOpt.beta .*tau).^Params.PulseOpt.n);
A_t((t < 0 | t>Trf)) = 0;
% disp( ['Average B1 of the pulse is:', num2str(mean(A_t))]) 

% Scaling Factor 
lambda = ((Params.PulseOpt.A0)^2 ./ (Params.PulseOpt.beta.*Params.PulseOpt.Q))^2;

% Frequency modulation function 
% Carrier frequency modulation function w(t):
omegaterm1 = sech((Params.PulseOpt.beta .* tau).^Params.PulseOpt.n);
omegaterm2 = cumtrapz(t,(omegaterm1.^2));
omega1 = -lambda*(omegaterm2 - omegaterm2(round(nSamples/2)))./(2*pi); % offset to allow for center at zero and rad/s to Hz
omega = lambda*omegaterm2;

% Phase modulation function phi(t):
phi = cumtrapz(tau, omega);

% Put together complex RF pulse waveform:
rf_pulse = A_t .* exp(1i .* phi);

%% Test 

% you need to find a way to integrate from t= 0 to each of the next time
% points
% This will give you nsamples-1 points. Use interp1 to interpolate this
% back to nSamples number of points. You will need to plot this as it
% involves extrapolation to make sure the limits at the ends make sense.
% Also try different interpolating functions. 

% omega = zeros(1, nSamples-1);
% for i = 1:nSamples-1
%     tau_range = t(1:i);
%     omega(i) = trapz(tau_range, (sech(Params.PulseOpt.beta * tau_range .^ Params.PulseOpt.n)).^2);
% end

%omega1 = interp1(tau(1:nSamples-1),omega, t, 'linear');
% omega1 = interp1(omega,t);

% Phase modulation function phi(t):
% phi = zeros(1, nSamples-1);
% for i = 1:nSamples-1
%     phi(i) = trapz(t, omega1);
% end
% phi = interp1(phi,t);



