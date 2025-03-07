function mt_sat_b1sim_CheckPulse(obj)

% Helper function to let you see what your pulse looks like.

alpha = obj.Prot.PulseSequenceParams.Mat(7);
Trf = obj.Prot.PulseSequenceParams.Mat(6)/1000;
shape = obj.options.SequenceSimulations_SatPulseShape;
delta = obj.Prot.PulseSequenceParams.Mat(2);  
PulseOpt.TBW = obj.Prot.PulseSequenceParams.Mat(8);
PulseOpt.bw = obj.Prot.PulseSequenceParams.Mat(9);

% Fix defaults set to 0. Let qMRLab decide:
if PulseOpt.TBW == 0
    PulseOpt.TBW = [];
end
if PulseOpt.bw == 0
    PulseOpt.bw = [];
end

Pulse = GetPulse(alpha, delta, Trf, shape, PulseOpt);

h = figure;
ViewPulse(Pulse,'b1')