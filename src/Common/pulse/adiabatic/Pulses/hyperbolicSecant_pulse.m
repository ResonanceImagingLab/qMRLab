function [rf_pulse, omega1, A_t, Params] = hyperbolicSecant_pulse(Trf, Params)

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
%   A(t) = A0 * sech(Beta*t)
%   omega1(t) = -mu*Beta*tanh(Beta*t)
%   A0 is the peak amplitude in microTesla
%   Beta is a frequency modulation parameter in rad/s
%   mu is a phase modulation parameter (dimensionless)
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
% Written by Christopher Rowley 2023


if ~exist('dispFigure','var') || isempty(dispFigure) || ~isfinite(dispFigure)
    dispFigure = 0;      
end


% Function to fill default values;
Params.PulseOpt = defaultHyperbolicSecParams(Params.PulseOpt);

nSamples = Params.PulseOpt.nSamples;  
t = linspace(0, Trf, nSamples);

% Amplitude
A_t =  Params.PulseOpt.A0* sech(Params.PulseOpt.beta* ( (t - Trf /2)).^Params.PulseOpt.n);
A_t((t < 0 | t>Trf)) = 0;
% disp( ['Average B1 of the pulse is:', num2str(mean(A_t))]) 


% Frequency modulation function 
% Carrier frequency modulation function w(t):
omega1 = -Params.PulseOpt.mu.*Params.PulseOpt.beta .* ...
            tanh(Params.PulseOpt.beta .* (t - Trf/2))./(2*pi); % 2pi to convert from rad/s to Hz

% Phase modulation function phi(t):
phi = Params.PulseOpt.mu .* log(sech(Params.PulseOpt.beta .* (t - Trf/2)) );

% Put together complex RF pulse waveform:
rf_pulse = A_t .* exp(1i .* phi);










































