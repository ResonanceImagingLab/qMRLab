function plotAdiabaticPulse(ax1,ax2, t, inv_pulse, Params)


%% Adiabatic Pulse Plot function 
%
% Called in adiabaticExample.m to check the amplitude and frequency
% modulation of your pulse 
%
%   inv_pulse denotes the adiabatic inversion pulse and can be found from
%   the individual pulse functions listed in getAdiabaticPulse.m 
%
% Written by Amie Demmans 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 20; 
%figure; 

% Plot Amplitude Function
    %ax1 = axes(adiabaticPulse);
    plot(ax1, t*1000, abs(inv_pulse), 'LineWidth', 3); 
    xlabel(ax1, 'Time(ms)'); 
    ylabel(ax1, 'B_1 (μT)')
    title(ax1,'Amplitude Function','FontWeight','normal');
    ax1.FontSize = fs;
    %set(gca, 'FontSize', fs);

% Plot Frequency Modulation Function
    %ax2 = axes(adiabaticPulse);
    %subplot(adiabaticPulse, 1,2,2); 
    plot(ax2, t*1000, imag(inv_pulse), 'LineWidth', 3);
    xlabel(ax2, 'Time(ms)'); 
    ylabel(ax2,'Frequency (Hz)');
    title(ax2,'Frequency Modulation function','FontWeight','normal');
    ax2.FontSize = fs;
    %set(gca, 'FontSize', fs);

% Scale display size     
%set(gcf,'Position',[100 100 1200 500])

% Set title of pulse 
%sgtitle(ax1, ax2, Params.Inv.shape, 'FontSize', fs+4,'FontWeight','bold')