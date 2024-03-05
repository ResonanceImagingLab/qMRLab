function PulseOpt = defaultLorentzParams(PulseOpt)

% Function designed to be used with the adiabatic pulse code
% Fills in default values if they are not user-defined

if(~isfield(PulseOpt,'beta') || isempty(PulseOpt.beta) || ~isfinite(PulseOpt.beta))
    % Default beta value in rad/s
    PulseOpt.beta = 99;       
end


if(~isfield(PulseOpt,'nSamples') || isempty(PulseOpt.nSamples) || ~isfinite(PulseOpt.nSamples))
    % Peak B1 of the pulse in microTesla
    PulseOpt.nSamples = 512;       
end

return;