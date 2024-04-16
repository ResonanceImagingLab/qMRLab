function AI_GUI
%% Adiabatic Inversion GUI 

% Open app figure
adibaticInv_GUI = uifigure("Name","Adiabatic Inversion Pulses");

g = uigridlayout(adibaticInv_GUI);
g.RowHeight = {22,22,'1x'};
g.ColumnWidth = {150,'1x'};
% Adding dropdown with pulse names 
AI_dropdown = uidropdown(adibaticInv_GUI,"Items",["Hs1","Lorentz","Gaussian","Hanning","Hsn","sin40"]);
