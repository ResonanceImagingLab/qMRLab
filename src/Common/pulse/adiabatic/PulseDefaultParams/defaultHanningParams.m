function PulseOpt = defaultHanningParams(PulseOpt)

% Function designed to be used with the adiabatic pulse code
% Fills in default values if they are not user-defined


if(~isfield(PulseOpt,'A0') || isempty(PulseOpt.A0) || ~isfinite(PulseOpt.A0))
    % Peak B1 of the pulse in microTesla
    PulseOpt.A0 = 14;       
end

if(~isfield(PulseOpt,'nSamples') || isempty(PulseOpt.nSamples) || ~isfinite(PulseOpt.nSamples))
    % Peak B1 of the pulse in microTesla
    PulseOpt.nSamples = 512;       
end


return;