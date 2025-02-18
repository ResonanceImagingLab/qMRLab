function Params = mt_sat_b1sim_getSeqParams(Params)

Params.lineshape = 'SuperLorentzian';
Params = mt_sat_b1sim_calcImagingParams(Params);

Params.WExcDur = 0.1/1000; % duration of water pulse
Params.PerfectSpoiling = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Params.numExcitation = 1; % forced to be 1. If doing something fancier, use ihMT module
Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'

% MT parameters that will be consistent:
Params.PulseOpt.bw = 0.3./Params.pulseDur; % override default Hann pulse shape.
Params.TD_MT =  0 ;   
Params.N_spin = 1;

if ~isfield(Params,'MTC') % if not defined, no MT
    Params.MTC = 1; % assume true if running the module
end

if Params.MTC
    % Allow different spoiling for MT
    if Params.GradientSpoiling 
        if ~isfield(Params,'GradientSpoilingStrength_MT')
            Params.GradientSpoilingStrength_MT = 20; % mT/m
        end
    end

    if ~isfield(Params,'A_g') || ~isfield(Params,'maxDephase') || ~isfield(Params,'G_t')
        % Calculate Spoiling moment from gradient
        Params = ihMT_calculateGradientSpoilingMoment( Params, 1 );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







