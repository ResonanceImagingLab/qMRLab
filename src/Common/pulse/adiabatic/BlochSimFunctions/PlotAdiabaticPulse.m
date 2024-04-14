function PlotAdiabaticPulse(t, inv_pulse, Params)

fs = 20; 
figure; 
tiledlayout(1,2)

% Plot Amplitude Function
    nexttile;
    plot(t*1000, abs(inv_pulse), 'LineWidth', 3); 
    xlabel('Time(ms)'); 
    ylabel('B_1 (μT)')
    title('Amplitude Function','FontWeight','normal');
    ax = gca; 
    ax.FontSize = fs;

% Plot Frequency Modulation Function
    
    nexttile; 
    plot(t*1000, imag(inv_pulse), 'LineWidth', 3);
    xlabel('Time(ms)'); 
    ylabel('Frequency (Hz)');
    title('Frequency Modulation function','FontWeight','normal');
    ax = gca; 
    ax.FontSize = fs;

% Scale display size     
set(gcf,'Position',[100 100 1200 500])

% Set title of pulse 
sgtitle(Params.Inv.shape, 'FontSize', fs+4,'FontWeight','bold')