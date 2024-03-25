function [rf_pulse, omega1, A_t, Params] = Lorentz_pulse( Trf, Params)

%   Lorentz_pulse Adiabatic Lorentz RF pulse function.
%   pulse = Lorentz_pulse(Trf, PulseOpt)
%
%   B1(t) = A(t) * exp( -1i *integral(omega1(t')) dt' )
%   where A(t) is the envelope, omega1 is the frequency sweep
%
%   Phase modulation is found from taking the integral of omega1(t)
%   Frequency modulation is time derivative of phi(t)
%
%   For the case of a Lorentz pulse:
%   A(t) = A_0 / (1 + Beta^2*t^2)
%   lambda = A_0^2/(Beta*Q)
%   omega1(t) = lambda(arctan(Beta*t)+Beta*t/(1+Beta^2*t^2))/2
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
%
%              de Graaf, R. A., & Nicolay, K. (1998). "Adiabatic water 
%              suppression using frequency selective excitation." Magnetic 
%              resonance in medicine, 40(5), 690-696.
%                   --> tau = (2*t/Trf)-1
%
%              Kupce, E. and Freeman, R (1995). "Optimized Adiabatic Pulses
%              for Wideband Spin Inversion." Journal of Magnetic Resonance
%              Imaging, Series A 118(2): 299-303.
%                  --> A_t, omega1, lambda equation 
%                  --> Q value of 4 
%
%              Tannus, A. Garwood, M. (1996). "Improved Performance of 
%              Frequency Swept Pulses Using Offset-Independent
%              Adiabaticity" Journal of Magnetic Resonance, 120(1),
%              133-137.
%                   --> Fig 1a and 1b. Show how width of amplitude and
%                   frequency vary with each pulse 
%                   --> A0 set to 20 as Lorentz pulse has the highest
%                   amplitude
%
% To be used with qMRlab
% Written by Christopher Rowley 2023 & Amie Demmans 2024


if ~exist('dispFigure','var') || isempty(dispFigure) || ~isfinite(dispFigure)
    dispFigure = 0;      
end

% Trf = 10;
% Params = Params.Inv;

% Function to fill default values;
Params.PulseOpt = defaultLorentzParams(Params.PulseOpt);

nSamples = Params.PulseOpt.nSamples;  
t = linspace(0, Trf, nSamples);

tau = ((2*t/Trf)-1);
%tau = (t - Trf /2);

% Amplitude
%A_t = Params.PulseOpt.A0./(1+(Params.PulseOpt.beta).*(tau.^2));
A_t = Params.PulseOpt.A0./(1+(Params.PulseOpt.beta.^2).*(tau.^2));
A_t((t < 0 | t>Trf)) = 0;
% disp( ['Average B1 of the pulse is:', num2str(mean(A_t))]) 

% Scaling Factor 
lambda = (Params.PulseOpt.A0)^2 ./ (Params.PulseOpt.beta.*Params.PulseOpt.Q);

% Frequency modulation function 
% Carrier frequency modulation function w(t)
omega1 = -lambda*(atan(Params.PulseOpt.beta.*tau)+(Params.PulseOpt.beta.*tau) ...
            /(1+(Params.PulseOpt.beta.^2).*(tau.^2)))./(2*pi); % convert rad/s to Hz

% Phase modulation function phi(t):
phi = lambda.*tau.*atan(Params.PulseOpt.beta.*tau)./2;

% Put together complex RF pulse waveform:
rf_pulse = A_t .* exp(1i .* phi);



%% Can do Bloch Sim to get inversion profile and display figure if interested:

% Params.NumPools = 1;
% BlochSimCallFunction(Params, rf_pulse, t, A_t, omega1);

% 
%     M_start = [0, 0, 0, 0, Params.M0a, Params.M0b]';
%     b1Rel = 0.5:0.1:1.5;
%     freqOff = -2000:200:2000;
%     [b1m, freqm] = ndgrid(b1Rel, freqOff);
% 
%     Mza = zeros(size(b1m));
%     Mzb = zeros(size(b1m));
% 
%     for i = 1:length(b1Rel)
%         for j = 1:length(freqOff)
% 
%             M_return = blochSimAdiabaticPulse( b1Rel(i)*rf_pulse, Params.Inv,  ...
%                             freqOff(j), Params, M_start, []);
% 
%             Mza(i,j) = M_return(5);
%             Mzb(i,j) = M_return(6);
%         end
%     end
% 
%     figure ('Name', 'Lorentz', 'NumberTitle', 'off'); 
%     tiledlayout(2,2)
%     nexttile; plot(t*1000, A_t, 'LineWidth', 3); 
%     xlabel('Time(ms)'); ylabel('B_1 (μT)')
%     title('Amplitude Function');ax = gca; ax.FontSize = 20;
% 
%     nexttile; plot(t*1000, omega1, 'LineWidth', 3);
%     xlabel('Time(ms)'); ylabel('Frequency (Hz)');
%     title('Frequency Modulation function');ax = gca; ax.FontSize = 20;
% 
%     nexttile; surf(b1m, freqm, Mza);
%     xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{za}');ax = gca; ax.FontSize = 20;
% 
%     nexttile; surf(b1m, freqm, Mzb);
%     xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{zb}');ax = gca; ax.FontSize = 20;
% 
%     set(gcf,'Position',[100 100 1200 1000])
% end

