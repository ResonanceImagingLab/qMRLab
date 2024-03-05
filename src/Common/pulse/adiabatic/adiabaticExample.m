%% Example  script Hsn
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);

Params.Inv.PulseOpt.beta = 672;  
Params.Inv.PulseOpt.n = 1; 
Params.Inv.PulseOpt.mu = 5;
Params.Inv.PulseOpt.A0 = 13.726;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10.24/1000;
Params.Inv.shape = 'hsn';
%Params.Inv.shape = 'lorentz';
Params.NumPools = 1; % or 2

%% Inversion
Params.dispFigure = 0;
[inv_pulse, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, ...
                                    Params.Inv);

% Check plot -> function for only plotting pulse (1 pool)
t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples);
    figure; tiledlayout(1,2)
    nexttile; plot(t*1000, abs(inv_pulse), 'LineWidth', 3); 
    xlabel('Time(ms)'); ylabel('B_1 (μT)')
    title('Amplitude Function');ax = gca; ax.FontSize = 20;
    
    nexttile; plot(t*1000, imag(inv_pulse), 'LineWidth', 3);
    xlabel('Time(ms)'); ylabel('Frequency (Hz)');
    title('Frequency Modulation function');ax = gca; ax.FontSize = 20;
set(gcf,'Position',[100 100 1200 500])

% bloch sim and return magnetization
BlochSimCallFunction

M_t = blochSimAdiabaticPulse( inv_pulse, Params.Inv ,...
                                  0, Params, M_t, []);
