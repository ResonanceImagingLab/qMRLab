function [Params, outputSamplingTable] = ihMT_getSeqParams(Params)

Params.lineshape = 'SuperLorentzian';
Params = ihMT_calcImagingParams(Params);

Params.WExcDur = 0.1/1000; % duration of water pulse
Params.PerfectSpoiling = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'

% MT parameters that will be consistent:
Params.PulseOpt.bw = 0.3./Params.pulseDur; % override default Hann pulse shape.
Params.TD_MT =  Params.TR_MT - Params.numSatPulse* (Params.pulseDur + Params.pulseGapDur) ;   
Params = ihMT_calcVariableImagingParams(Params);


% Other image parameters
Params.NumLines = 216;
Params.NumPartitions = 192; 
Params.Slices = 176;
Params.Grappa = 1;
Params.ReferenceLines = 32;
Params.AccelerationFactor = 2;
Params.Segments = []; 
Params.TurboFactor = Params.numExcitation- Params.DummyEcho;
Params.ellipMask = 1;

[outputSamplingTable, ~, Params.Segments] = ihMT_Step1_calculateKspaceSampling_v3 (Params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







