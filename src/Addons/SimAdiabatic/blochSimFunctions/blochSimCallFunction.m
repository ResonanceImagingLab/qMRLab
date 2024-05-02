function blochSimCallFunction(rf_pulse, Params)
%% Can do Bloch Sim to get inversion profile and display figure if interested:

fs = 18;

if Params.NumPools == 1
    M_start = [0, 0, Params.M0a]';
    b1Rel = 0.5:0.1:1.5; %B1 inhomogeneity
    freqOff = -2000:200:2000; %B0 inhomogeneity
    [b1m, freqm] = ndgrid(b1Rel, freqOff); % generates grid 
    
    Mza = zeros(size(b1m));
    
    for i = 1:length(b1Rel)
        for j = 1:length(freqOff)
        
            M_return = blochSimAdiabaticPulse_1pool( b1Rel(i)*rf_pulse, Params.Inv,  ...
                            freqOff(j), Params, M_start, []);

            Mza(i,j) = M_return(3);
        end
    end

     figure; tiledlayout(1,1)
    % nexttile; plot(t*1000, A_t, 'LineWidth', 3); 
    % xlabel('Time(ms)'); ylabel('B_1 (μT)')
    % title('Amplitude Function','FontWeight','normal');ax = gca; ax.FontSize = fs;
    
    % nexttile; plot(t*1000, omega1, 'LineWidth', 3);
    % xlabel('Time(ms)'); ylabel('Frequency (Hz)');
    % title('Frequency Modulation function','FontWeight','normal');ax = gca; ax.FontSize = fs;
    
    nexttile; surf(b1m, freqm, Mza,'FaceColor','interp'); 
    xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{za}');ax = gca; ax.FontSize = fs;
    view(-55, 10)
    
    set(gcf,'Position',[100 100 700 500])

    sgtitle(Params.Inv.shape, 'FontSize', fs+4,'FontWeight','bold')

elseif Params.NumPools == 2
    M_start = [0, 0, 0, 0, Params.M0a, Params.M0b]';
    b1Rel = 0.5:0.1:1.5;
    freqOff = -2000:200:2000;
    [b1m, freqm] = ndgrid(b1Rel, freqOff);
    
    Mza = zeros(size(b1m));
    Mzb = zeros(size(b1m));
    
    for i = 1:length(b1Rel)
        for j = 1:length(freqOff)
        
            M_return = blochSimAdiabaticPulse_2pool( b1Rel(i)*rf_pulse, Params.Inv,  ...
                            freqOff(j), Params, M_start, []);

            Mza(i,j) = M_return(5);
            Mzb(i,j) = M_return(6);
        end
    end

    figure; tiledlayout(1,2)
    % nexttile; plot(t*1000, A_t, 'LineWidth', 3); 
    % xlabel('Time(ms)'); ylabel('B_1 (μT)')
    % title('Amplitude Function','FontWeight','normal');ax = gca; ax.FontSize = fs;
    % 
    % nexttile; plot(t*1000, omega1, 'LineWidth', 3);
    % xlabel('Time(ms)'); ylabel('Frequency (Hz)');
    % title('Frequency Modulation function','FontWeight','normal');ax = gca; ax.FontSize = fs;
    
    nexttile; surf(b1m, freqm, Mza, 'FaceColor','interp');
    xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{za}');ax = gca; ax.FontSize = fs;
    
    nexttile; surf(b1m, freqm, Mzb,'FaceColor','interp');
    xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{zb}');ax = gca; ax.FontSize = fs;
    
    set(gcf,'Position',[100 100 1000 500])

    sgtitle(Params.Inv.shape, 'FontSize', fs+4,'FontWeight','bold')

else
    error('Define Params.NumPools to be = 1 or 2;')
end

return; 
