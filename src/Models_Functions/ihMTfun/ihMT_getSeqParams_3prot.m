function [Params, outputSamplingTable] = ihMT_getSeqParams_3prot(Params)

%Params.B0 = 3;
%Params.MTC = 1; % Magnetization Transfer Contrast
%Params.MTC = obj.Prot.PulseSequenceParams.Mat(1);
%Params.TissueType = 'GM';
%Params = ihMT_defaultCortexTissueParams(Params);
Params.lineshape = 'SuperLorentzian';
Params = ihMT_calcImagingParams(Params);

Params.WExcDur = 0.1/1000; % duration of water pulse
Params.echospacing = 7.66/1000;
Params.PerfectSpoiling = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Params listed in this section need to change into this format:

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Params.TurboFactor = obj.Prot.PulseSequenceParams.Mat(6);
Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'

% MT parameters that will be consistent:
Params.SatPulseShape = 'gausshann';
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
% Params for TF = 80 in the two diff ways 

% Params.delta = 8000;
% Params.flipAngle = 7; % excitation flip angle water.
% Params.TR = 1.14; % total repetition time = MT pulse train and readout.
% Params.numSatPulse = 6;
% Params.TurboFactor = 80;
% Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
% Params.satFlipAngle = 11.6; % microTesla  
% Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
% Params.DummyEcho = 2;
% Params.boosted = 1;
% Params.satTrainPerBoost = 9;
% Params.TR_MT = 0.06; 


%Params.delta = obj.Prot.PulseSequenceParams.Mat(2);  % delta 
%Params.flipAngle = obj.Prot.PulseSequenceParams.Mat(3);
%Params.TR = obj.Prot.PulseSequenceParams.Mat(4);
%Params.numSatPulse = obj.Prot.PulseSequenceParams.Mat(5);
%Params.TurboFactor = obj.Prot.PulseSequenceParams.Mat(6);
%Params.pulseDur = obj.Prot.PulseSequenceParams.Mat(7);
%Params.satFlipAngle = obj.Prot.PulseSequenceParams.Mat(8);
%Params.pulseGapDur = obj.Prot.PulseSequenceParams.Mat(9);
%Params.DummtEcho = obj.Prot.PulseSequenceParams.Mat(10);
%Params.boosted = obj.Prot.PulseSequenceParams.Mat(11);
%Params.satTrainPerBoost = obj.Prot.PulseSequenceParams.Mat(12);
%Params.TR_MT = obj.Prot.PulseSequenceParams.Mat(13);








