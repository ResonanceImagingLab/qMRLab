%% Adiabatic Inversion (AI) GUI 
% Using dot notation as that is what is followed in GUI development Matlab 

% Open app figure
adiabaticInv_GUI = uifigure("Name","Adiabatic Inversion Pulses");
adiabaticInv_GUI.Position = [100 100 1000 600];

% Use grid layout manager to manage position and size of UI components 
g = uigridlayout(adiabaticInv_GUI);
g.RowHeight = {22,22,'1x'}; % insert later when more components added 
g.ColumnWidth = {150,'1x'}; % insert later when more components added 

% Create label for pulse shape drop down
pulse_lbl = uilabel(g);
pulse_lbl.Text = "Pulse:";
pulse_lbl.Layout.Row = 1;
pulse_lbl.Layout.Column = 1;

% Adding dropdown with pulse names 
AI_dropdown = uidropdown(g);
AI_dropdown.Items = ["Hs1","Lorentz","Gauss","Hanning","Hsn","Sin40"];
AI_dropdown.Layout.Row = 2;
AI_dropdown.Layout.Column = 1;
AI_dropdown.Value = "Hs1";

% Create label for B0 magnet size dropdown 
AI_magnet_lbl = uilabel(g);
AI_magnet_lbl.Text = "Magnet (Tesla):";
AI_magnet_lbl.Layout.Row = 1; 
AI_magnet_lbl.Layout.Column = 2;

% Create dropdown for B0 magnet size 
AI_magnet_dropdown = uidropdown(g);
AI_magnet_dropdown.Items = ["1.5","3.0","7.0"];
AI_magnet_dropdown.Value = "3.0";
AI_magnet_dropdown.Layout.Row = 2;
AI_magnet_dropdown.Layout.Column = [2,3]; 

% Create labels for edit fields 
% Beta
Beta_lbl = uilabel(g);
Beta_lbl.Text = "Beta (rad/s):";
Beta_lbl.Layout.Row = 3;
Beta_lbl.Layout.Column = 1;
% n 
n_lbl = uilabel(g);
n_lbl.Text = "n:";
n_lbl.Layout.Row = 5;
n_lbl.Layout.Column = 1; 
%A0 
A0_lbl = uilabel(g);
A0_lbl.Text = "A0 (microTesla):";
A0_lbl.Layout.Row = 7;
A0_lbl.Layout.Column = 1; 
% nSamples 
nSamples_lbl = uilabel(g);
nSamples_lbl.Text = "nSamples:";
nSamples_lbl.Layout.Row = 9;
nSamples_lbl.Layout.Column = 1; 
% Q
Q_lbl = uilabel(g);
Q_lbl.Text = "Q:";
Q_lbl.Layout.Row = 11;
Q_lbl.Layout.Column = 1;
% Trf 
Trf_lbl = uilabel(g);
Trf_lbl.Text = 'Trf (s):';
Trf_lbl.Layout.Row = 13;
Trf_lbl.Layout.Column = 1;

% Set edit fields for each parameter 
% Beta 
edit_Beta = uieditfield(g, 'numeric');
edit_Beta.Tag = 'beta';
edit_Beta.Layout.Row = 4;
edit_Beta.Layout.Column = 1; 
% n
edit_n = uieditfield(g,'numeric');
edit_n.Tag = 'n';
edit_n.Layout.Row = 6; 
edit_n.Layout.Column = 1; 
% A0 
edit_A0 = uieditfield(g,'numeric');
edit_A0.Tag = 'A0';
edit_A0.Layout.Row = 8; 
edit_A0.Layout.Column = 1;
% nSamples 
edit_nSamples = uieditfield(g,'numeric');
edit_nSamples.Tag = 'nSamples';
edit_nSamples.Layout.Row = 10; 
edit_nSamples.Layout.Column = 1;
% Q
edit_Q = uieditfield(g,'numeric');
edit_Q.Tag = 'Q';
edit_Q.Layout.Row = 12; 
edit_Q.Layout.Column = 1; 
% Trf 
edit_Trf = uieditfield(g,'numeric');
edit_Trf.Tag = 'Trf';
edit_Trf.Layout.Row = 14;
edit_Trf.Layout.Column = 1;

% Adding uiaxes for plotting the functions 
% --> May need to set 4 of these for the 4 plots in a 2 pool case
% Later see if there is a way so it is not having to set individually 
ax1 = uiaxes(g);
ax1.Layout.Row = [3,14];
ax1.Layout.Column = 2;

ax2 = uiaxes(g);
ax2.Layout.Row = [3,14];
ax2.Layout.Column = 3;

ax1.Parent = g;
ax2.Parent = g;


%setDefaultParams(g, AI_magnet_dropdown, AI_dropdown, edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2);

% Set callback for dropdown and edit fields 
AI_dropdown.ValueChangedFcn = @(dropdown, event) dropdownMenu_Callback(dropdown, event, g, AI_magnet_dropdown, edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2);
edit_Beta.ValueChangedFcn = @(editField, event) editFieldChangedCallback(editField, event, g, AI_magnet_dropdown, edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2);
edit_n.ValueChangedFcn = @(editField, event) editFieldChangedCallback(editField, event, g, AI_magnet_dropdown, edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2);
edit_A0.ValueChangedFcn = @(editField, event) editFieldChangedCallback(editField, event, g, AI_magnet_dropdown, edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2);
edit_nSamples.ValueChangedFcn = @(editField, event) editFieldChangedCallback(editField, event, g, AI_magnet_dropdown, edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2);
edit_Q.ValueChangedFcn = @(editField, event) editFieldChangedCallback(editField, event, g, AI_magnet_dropdown, edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2);
edit_Trf.ValueChangedFcn = @(editField, event) editFieldChangedCallback(editField, event, g, AI_magnet_dropdown, edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2);

% Function to set default parameters for selected pulse 
function dropdownMenu_Callback(dropdown, event, g, AI_magnet_dropdown, edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2)

%pulse = AI_dropdown.Value;
pulse = dropdown.Value;
%magnetSize = str2double(AI_magnet_dropdown.Value);

%Params = struct();
%Params.B0 = magnetSize;
Params.B0 = str2double(AI_magnet_dropdown.Value);
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
                Params = struct(Params);
         end

 % Setting values of edit fields 
    setEditField(edit_Beta, Params.beta);
    setEditField(edit_n, Params.n);
    setEditField(edit_A0, Params.A0);
    setEditField(edit_nSamples,Params.nSamples);
    setEditField(edit_Q, Params.Q);
    setEditField(edit_Trf, Params.Trf);

    % Params.beta = edit_Beta.Value;
    % Params.n = edit_n.Value;
    % Params.A0 = edit_A0.Value; 
    % Params.nSamples = edit_nSamples.Value;
    % Params.Q = edit_Q.Value;
    % Params.Trf = edit_Trf.Value; 
 % Plotting Functions 
        [inv_pulse, omega1, A_t, Params] = getAdiabaticPulse( Params.Trf, pulse, Params);
        t = linspace(0, Params.Trf, Params.nSamples);
        plotAdiabaticPulse(ax1, ax2, t, inv_pulse,Params);
end

% Function to allow edit field changing to apply to plots 
function editFieldChangedCallback(editField, event, AI_dropdown, AI_magnet_dropdown,edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2)

dropdownMenu_Callback(AI_dropdown, event, AI_magnet_dropdown, edit_Beta, edit_n,edit_A0, edit_nSamples, edit_Q, edit_Trf, ax1, ax2);
end


function setEditField(editField, value)
    % Set the value of the edit field
    editField.Value = value;
   end

















