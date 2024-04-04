function PulseOpt = defaultHsnParams(PulseOpt)

% Function designed to be used with the adiabatic pulse code
% Fills in default values if they are not user-defined
% If n is changed from 8 these params will also need to change 

if(~isfield(PulseOpt,'beta') || isempty(PulseOpt.beta) || ~isfinite(PulseOpt.beta))
    % Default beta value in rad/s
    PulseOpt.beta = 250;       
end

if(~isfield(PulseOpt,'n') || isempty(PulseOpt.n) || ~isfinite(PulseOpt.n))
    % Default sech exponent
    PulseOpt.n = 8;       
end

if(~isfield(PulseOpt,'A0') || isempty(PulseOpt.A0) || ~isfinite(PulseOpt.A0))
    % Peak B1 of the pulse in microTesla
    PulseOpt.A0 = 17.5;       
end

if(~isfield(PulseOpt,'nSamples') || isempty(PulseOpt.nSamples) || ~isfinite(PulseOpt.nSamples))
    % Peak B1 of the pulse in microTesla
    PulseOpt.nSamples = 512;       
end

if(~isfield(PulseOpt,'Q') || isempty(PulseOpt.Q) || ~isfinite(PulseOpt.Q))
    % Adiabaticity Factor
    PulseOpt.nSamples = 4e-4;       
end

return;