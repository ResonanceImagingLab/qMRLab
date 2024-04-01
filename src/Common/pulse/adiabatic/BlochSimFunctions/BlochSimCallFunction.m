function BlochSimCallFunction(rf_pulse, t, A_t, omega1,Params)
%% Can do Bloch Sim to get inversion profile and display figure if interested:



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

    figure; tiledlayout(1,3)
    nexttile; plot(t*1000, A_t, 'LineWidth', 3); 
    xlabel('Time(ms)'); ylabel('B_1 (μT)')
    title('Amplitude Function');ax = gca; ax.FontSize = 20;
    
    nexttile; plot(t*1000, omega1, 'LineWidth', 3);
    xlabel('Time(ms)'); ylabel('Frequency (Hz)');
    title('Frequency Modulation function');ax = gca; ax.FontSize = 20;
    
    nexttile; surf(b1m, freqm, Mza,'FaceColor','interp'); 
    xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{za}');ax = gca; ax.FontSize = 20;
    view(-55, 10)
    
    set(gcf,'Position',[100 100 1200 500])

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

    figure; tiledlayout(2,2)
    nexttile; plot(t*1000, A_t, 'LineWidth', 3); 
    xlabel('Time(ms)'); ylabel('B_1 (μT)')
    title('Amplitude Function');ax = gca; ax.FontSize = 20;
    
    nexttile; plot(t*1000, omega1, 'LineWidth', 3);
    xlabel('Time(ms)'); ylabel('Frequency (Hz)');
    title('Frequency Modulation function');ax = gca; ax.FontSize = 20;
    
    nexttile; surf(b1m, freqm, Mza, 'FaceColor','interp');
    xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{za}');ax = gca; ax.FontSize = 20;
    
    nexttile; surf(b1m, freqm, Mzb,'FaceColor','interp');
    xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{zb}');ax = gca; ax.FontSize = 20;
    
    set(gcf,'Position',[100 100 1200 1000])

else
    error('Define Params.NumPools to be = 1 or 2;')
end

return; 
