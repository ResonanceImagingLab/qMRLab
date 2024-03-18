function PulseOpt = defaultGaussParams(PulseOpt)

% Function designed to be used with the adiabatic pulse code
% Fills in default values if they are not user-defined

if(~isfield(PulseOpt,'beta') || isempty(PulseOpt.beta) || ~isfinite(PulseOpt.beta))
    % Default beta value in rad/s
    PulseOpt.beta = 20;       
end

if(~isfield(PulseOpt,'A0') || isempty(PulseOpt.A0) || ~isfinite(PulseOpt.A0))
    % Peak B1 of the pulse in microTesla
    PulseOpt.A0 = 17.5;       
end

if(~isfield(PulseOpt,'nSamples') || isempty(PulseOpt.nSamples) || ~isfinite(PulseOpt.nSamples))
    % Default number of samples taken based on machine properties
    PulseOpt.nSamples = 512;       
end

if(~isfield(PulseOpt,'Q') || isempty(PulseOpt.Q) || ~isfinite(PulseOpt.Q))
    % Adiabaticity factor 
    PulseOpt.Q = 4;       
end

return;