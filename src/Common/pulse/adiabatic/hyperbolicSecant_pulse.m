function [A_t, PulseOpt] = hyperbolicSecant_pulse(t, Trf, PulseOpt)

%   hyperbolicSecant_pulse Adiabatic hyperbolic secant RF pulse function.
%   pulse = hyperbolicSecant_pulse(t, Trf, PulseOpt)
%
%   B1(t) = A(t) * exp( -1i *integral(omega1(t')) dt' )
%   where A(t) is the envelope, omega1 is the frequency sweep
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
% 
%   Reference: Matt A. Bernstein, Kevin F. Kink and Xiaohong Joe Zhou.
%              Handbook of MRI Pulse Sequences, pp. 110, Eq. 4.10, (2004)
%
%              Tannús, A. and M. Garwood (1997). "Adiabatic pulses." 
%              NMR in Biomedicine 10(8): 423-434.
%
%   See also GETPULSE, VIEWPULSE.
%
% To be used with qMRlab
% Written by Christopher Rowley 2023

if (nargin < 3); PulseOpt = struct; end

% Function to fill default values;
PulseOpt = defaultHyperbolicSecParams(PulseOpt);


% PulseOpt.A0 = 12;
nSamples = pulseOpt.nSamples;  
t = linspace(0, Trf, nSamples);

% Amplitude
A_t =  PulseOpt.A0* sech(PulseOpt.beta* ( (t - Trf/2)).^PulseOpt.n);
A_t((t < 0 | t>Trf)) = 0;
disp( ['Average B1 of the pulse is:', num2str(mean(A_t))]) 


% Frequency modmodulation function 
% Carrier frequency modulation function w(t):
omega1 = -PulseOpt.mu.*PulseOpt.beta .* ...
            tanh(PulseOpt.beta .* (t - Trf/2))./(2*pi); % 2pi to convert from rad/s to Hz

% Phase modulation function phi(t):
phi = PulseOpt.mu .* log(sech(PulseOpt.beta .* (t - Trf/2)) );

% Put together complex RF pulse waveform:
rf_pulse = A_t .* exp(1i .* phi);


%% Bloch Sim to get inversion profile
b1Rel = linspace(0.5, 1.5, 10);
freqOff = -2000:200:2000;
[b1m, freqm] = ndgrid(b1Rel, freqOff);

Mza = zeros(size(b1m));
Mzb = zeros(size(b1m));

for i = 1:length(b1Rel)
    for j = 1:length(freqOff)

        M_return = blochSimAdiabaticPulse( b1Rel(i)* rf_pulse,...
            Trf, PulseOpt, freqOff(j));

        Mza(i,j) = M_return(5);
        Mzb(i,j) = M_return(6);
    end
end

figure; tiledlayout(2,2)
nexttile; plot(t, A_t, 'LineWidth', 3); 
title('Amplitude Function');ax = gca; ax.FontSize = 20;

nexttile; plot(t, omega1, 'LineWidth', 3);
title('Frequency Modulation function');ax = gca; ax.FontSize = 20;

nexttile; surf(b1m, freqm, Mza);
xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{za}');ax = gca; ax.FontSize = 20;

nexttile; surf(b1m, freqm, Mzb);
xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{zb}');ax = gca; ax.FontSize = 20;

set(gcf,'Position',[100 100 1200 1000])

return; 








































