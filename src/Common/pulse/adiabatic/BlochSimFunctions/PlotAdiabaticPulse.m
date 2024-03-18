function PlotAdiabaticPulse(t, inv_pulse)


figure; 
tiledlayout(1,2)

% Plot Amplitude Function
    nexttile;
    plot(t*1000, abs(inv_pulse), 'LineWidth', 3); 
    xlabel('Time(ms)'); 
    ylabel('B_1 (μT)')
    title('Amplitude Function');
    ax = gca; 
    ax.FontSize = 20;

% Plot Frequency Modulation Function
    
    nexttile; 
    plot(t*1000, imag(inv_pulse), 'LineWidth', 3);
    xlabel('Time(ms)'); 
    ylabel('Frequency (Hz)');
    title('Frequency Modulation function');
    ax = gca; 
    ax.FontSize = 20;

% Scale display size     
set(gcf,'Position',[100 100 1200 500])