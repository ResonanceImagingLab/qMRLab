%function AI_GUI
%% Adiabatic Inversion (AI) GUI 

% Open app figure
adibaticInv_GUI = uifigure("Name","Adiabatic Inversion Pulses");

% Use grid layout manager to manage position and size of UI components 
g = uigridlayout(adibaticInv_GUI);
g.RowHeight = {22,22,'1x'}; % insert later when more components added 
g.ColumnWidth = {150,'1x'}; % insert later when more components added 

% Create label for pulse shape drop down
pulse_lbl = uilabel(g);
pulse_lbl.Text = "Pulse:";
pulse_lbl.Layout.Row = 1;
pulse_lbl.Layout.Column = 1;

% Adding dropdown with pulse names 
% Using dot notation as that is what is followed in GUI development Matlab 
AI_dropdown = uidropdown(g);
AI_dropdown.Items = ["Hs1","Lorentz","Gauss","Hanning","Hsn","Sin40"];
AI_dropdown.Layout.Row = 2;
AI_dropdown.Layout.Column = 1;
%AI_dropdown.Value = "Hs1";


% Create labels for edit fields 
% Beta
Betalbl = uilabel(g);
Betalbl.Text = "Beta:";
Betalbl.Layout.Row = 3;
Betalbl.Layout.Column = 1;
% n 
nlbl = uilabel(g);
nlbl.Text = "n:";
nlbl.Layout.Row = 5;
nlbl.Layout.Column = 1; 
%A0 
A0lbl = uilabel(g);
A0lbl.Text = "A0:";
A0lbl.Layout.Row = 7;
A0lbl.Layout.Column = 1; 
% nSamples 
nSampleslbl = uilabel(g);
nSampleslbl.Text = "nSamples:";
nSampleslbl.Layout.Row = 9;
nSampleslbl.Layout.Column = 1; 
% Q
Qlbl = uilabel(g);
Qlbl.Text = "Q:";
Qlbl.Layout.Row = 11;
Qlbl.Layout.Column = 1;

% Set edit fields for each parameter 
% Beta 
editBeta = uieditfield(g, 'numeric');
editBeta.Tag = 'beta';
editBeta.Layout.Row = 4;
editBeta.Layout.Column = 1; 
% n
editn = uieditfield(g,'numeric');
editn.Tag = 'n';
editn.Layout.Row = 6; 
editn.Layout.Column = 1; 
% A0 
editA0 = uieditfield(g,'numeric');
editA0.Tag = 'A0';
editA0.Layout.Row = 8; 
editA0.Layout.Column = 1;
% nSamples 
editnSamples = uieditfield(g,'numeric');
editnSamples.Tag = 'nSamples';
editnSamples.Layout.Row = 10; 
editnSamples.Layout.Column = 1;
% Q
editQ = uieditfield(g,'numeric');
editQ.Tag = 'Q';
editQ.Layout.Row = 12; 
editQ.Layout.Column = 1; 

% Create label for magnet size dropdown 
AI_magnet_lbl = uilabel(g);
AI_magnet_lbl.Text = "Magnet (Tesla):";
AI_magnet_lbl.Layout.Row = 1; 
AI_magnet_lbl.Layout.Column = 2;

% Adding dropdown for Tesla magnet size 
AI_magnet = uidropdown(g);
AI_magnet.Items = ["3","1.5","7"];
AI_magnet.Layout.Row = 2;
AI_magnet.Layout.Column = 2;
AI_magnet.Value = '3';

% Set callback for dropdown menu
AI_dropdown.ValueChangedFcn = @(dropdown, event) dropdownMenu_Callback(dropdown, event, g, editBeta, editn,editA0, editnSamples, editQ);


% Adding dropdown for Tesla magnet size 
% AI_magnet = uidropdown(g);
% AI_magnet.Items = ["3.0","1.5","7T"];
% AI_magnet.Layout.Row = 1;
% AI_magnet.Layout.Column = 2;
% AI_magnet.Value = '3.0';


% Function to set default values based on dropdown selection  
function dropdownMenu_Callback(dropdown, event, g, editBeta, editn,editA0, editnSamples, editQ)
 
        pulse = dropdown.Value;
        
        Params.B0 = 3;
        Params.TissueType = 'WM'; % white matter 
        Params = AI_defaultTissueParams(Params);

         switch pulse 
            case "Hs1"
                Params = AI_defaultHs1Params(Params);
            case "Lorentz"
                Params = AI_defaultLorentzParams(Params);
            case "Gauss"
                Params = AI_defaultGaussParams(Params);
            case "Hanning"
                Params = AI_defaultHanningParams(Params);
            case "Hsn"
                Params = AI_defaultHsnParams(Params);
            case "Sin40"
                Params = AI_defaultSin40Params(Params);
            otherwise 
                Params= struct(Params);
         end

 % Setting values of edit fields 
    setEditField(editBeta, Params.beta);
    setEditField(editn, Params.n);
    setEditField(editA0, Params.A0);
    setEditField(editnSamples,Params.nSamples);
    setEditField(editQ, Params.Q);
    end

function setEditField(editField, value)
    % Set the value of the edit field
    editField.Value = value;
   end



















