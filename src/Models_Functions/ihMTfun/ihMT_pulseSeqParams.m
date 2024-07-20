function Params = ihMT_pulseSeqParams(Params)

if(~isfield(Params,'MTC') || isempty(Params.MTC) || ~isfinite(Params.MTC))
    % Default magnetization transfer contrast 
    Params.MTC = 1;       
end

if(~isfield(Params,'delta') || isempty(Params.delta) || ~isfinite(Params.delta))
    % Default ...
    Params.delta = 8000;       
end

if(~isfield(Params,'flipAngle') || isempty(Params.flipAngle) || ~isfinite(Params.flipAngle))
    % Default excitation flip angle of water 
    Params.flipAngle = 7;       
end

if(~isfield(Params,'TR') || isempty(Params.TR) || ~isfinite(Params.TR))
    % Default total repitition time = MT pulse train and readout 
    Params.TR = 1.14;       
end

if(~isfield(Params,'numSatPulse') || isempty(Params.numSatPulse) || ~isfinite(Params.numSatPulse))
    % Default beta value in rad/s (modulation angular frequency)
    Params.numSatPulse = 6;       
end

if(~isfield(Params,'TurboFactor') || isempty(Params.TurboFactor) || ~isfinite(Params.TurboFactor))
    % Default turbo factor 
    Params.TurboFactor = 80;       
end

if(~isfield(Params,'pulseDur') || isempty(Params.pulseDur) || ~isfinite(Params.pulseDur))
    % Duration of one mT pulse in seconds 
    Params.pulseDur = 0.768/1000;       
end

if(~isfield(Params,'satFlipAngle') || isempty(Params.satFlipAngle) || ~isfinite(Params.satFlipAngle))
    % Default saturation flip angle 
    Params.satFlipAngle = 11.6;       
end

if(~isfield(Params,'pulseGapDur') || isempty(Params.pulseGapDur) || ~isfinite(Params.pulseGapDur))
    % Default ms gap between MT pulses in train 
    Params.pulseGapDur = 0.3/1000;       
end

if(~isfield(Params,'DummyEcho') || isempty(Params.DummyEcho) || ~isfinite(Params.DummyEcho))
    % Default # of dummy echoes??
    Params.DummyEcho = 2;       
end

if(~isfield(Params,'boosted') || isempty(Params.boosted) || ~isfinite(Params.boosted))
    % Default low duty cycle
    Params.boosted = 1;       
end

if(~isfield(Params,'satTrainPerBoost') || isempty(Params.satTrainPerBoost) || ~isfinite(Params.satTrainPerBoost))
    % Default ...
    Params.satTrainPerBoost = 9;       
end

if(~isfield(Params,'TR_MT') || isempty(Params.TR_MT) || ~isfinite(Params.TR_MT))
    % Default ...
    Params.TR_MT = 0.06;       
end














