classdef vfa_t1_spoil < AbstractModel
% vfa_t1: Compute a T1 map using Variable Flip Angle
% This module is a port of the hMRI toolbox code to qMRlab
%
% Assumptions:
%
% Inputs:
%   VFAData         Spoiled Gradient echo data, 4D volume with different flip angles in time dimension
%   (B1map)         Normalized transmit excitation field map (B1+). B1+ is defined 
%                   as a  normalized multiplicative factor such that:
%                   FA_actual = B1+ * FA_nominal. (OPTIONAL).
%   (Mask)          Binary mask to accelerate the fitting. (OPTIONAL)
%
% Outputs:
%   T1              Longitudinal relaxation time [s]
%   M0              Equilibrium magnetization
%
% Protocol:
%   VFAData Array [nbFA x 2]:
%       [FA1 TR1; FA2 TR2;...]      flip angle [degrees] TR [s]
%
% Options:
%   None
%
% Example of command line usage:
%   Model = vfa_t1;  % Create class from model
%   Model.Prot.VFAData.Mat=[3 0.015; 20 0.015]; %Protocol: 2 different FAs
%   data = struct;  % Create data structure
%   data.VFAData = load_nii_data('VFAData.nii.gz');
%   data.B1map = load_nii_data('B1map.nii.gz');
%   FitResults = FitData(data,Model); %fit data
%   FitResultsSave_mat(FitResults);
%
%   For more examples: <a href="matlab: qMRusage(vfa_t1);">qMRusage(vfa_t1)</a>
%
%
% Author: Christopher ROwley, 2025
%
% References:
%   Please cite the following if you use this module:
%     Fram, E.K., Herfkens, R.J., Johnson, G.A., Glover, G.H., Karis, J.P.,
%     Shimakawa, A., Perkins, T.G., Pelc, N.J., 1987. Rapid calculation of
%     T1 using variable flip angle gradient refocused imaging. Magn. Reson.
%     Imaging 5, 201?208
%   In addition to citing the package:
%     Karakuzu A., Boudreau M., Duval T.,Boshkovski T., Leppert I.R., Cabana J.F., 
%     Gagnon I., Beliveau P., Pike G.B., Cohen-Adad J., Stikov N. (2020), qMRLab: 
%     Quantitative MRI analysis, under one umbrella doi: 10.21105/joss.02343
%
% Specific References for the hMRI toolbox:
% Tabelow, K., Balteau, E., Ashburner, J., Callaghan, M. F., Draganski, B.,
% Helms, G., Kherif, F., Leutritz, T., Lutti, A., Phillips, C., Reimer, E.,
% Ruthotto, L., Seif, M., Weiskopf, N., Ziegler, G., Mohammadi, S., 2019. 
% hMRI – A toolbox for quantitative MRI in neuroscience and clinical research. 
% Neuroimage 194, 191-210. https://doi.org/10.1016/j.neuroimage.2019.01.029

properties (Hidden=true)
 onlineData_url = 'https://osf.io/7wcvh/download?version=3';  
end

properties
    MRIinputs = {'LowFlipAngle','HighFlipAngle','B1map','Mask'};
    xnames = {'S0','T1'};
    voxelwise = 0;
    
    % Protocol
    Prot  = struct('VFAProtocol',struct('Format',{{'FlipAngle' 'TR (ms)'}},...
                     'Mat', [3 15; 20 15]), ...
                    'SpoilingCorrection', struct('Format',{{ ...
                    'RFspoilingIncrement(deg)'; 'SpoilingGradientDuration(ms)';...
                    'SpoilingGradientAmplitude(mT/m)'; 'SimT1_Lower(ms)'; ...
                    'SimT1_Upper(ms)'; 'AverageT2(ms)'; 'B1_Lower(rel)'; ...
                    'B1_Upper(rel)'; 'DiffusionCoefficient(um2/ms)'}}, ...
                    'Mat', [ 50; 3.38; 24; 500; 2500; 80; 0.5; 1.5; 0.8])); % You can define a default protocol here.

     buttons = {'PANEL', 'IncompleteSpoiling',4,...
    'smallAngleApprox',{'True','False'}, ...
    'ScannerVendor', {'Siemens', 'GE', 'Philips'}, ...
    'Run Spoiling Simulations','pushbutton',...
    'Load Spoiling Simulations Results','pushbutton'};


    % Model options
    options= struct(); % structure filled by the buttons. Leave empty in the code
    previousOptions = struct()
end

methods (Hidden=true)
% Hidden methods goes here.
end

methods

    function obj = vfa_t1_spoil()
        obj.options = button2opts(obj.buttons);
        obj.previousOptions = obj.options;
    end

    function checkfields = checkupdatedfields(obj)
        if (~isequal(obj.options.IncompleteSpoiling_ScannerVendor, obj.previousOptions.IncompleteSpoiling_ScannerVendor))
            % Flag to auto update the spoiling increment
            checkfields = 1; 
        elseif (~isequal(obj.options.IncompleteSpoiling_RunSpoilingSimulations, obj.previousOptions.IncompleteSpoiling_RunSpoilingSimulations)||...
                 ~isequal(obj.options.IncompleteSpoiling_LoadSpoilingSimulationsResults, obj.previousOptions.IncompleteSpoiling_LoadSpoilingSimulationsResults))
            % This is current set up to run M0B_R1obs, but need to run
            % the spoling
            checkfields = 2;

        else
            checkfields = 0;
        end 
    end 

    function obj = UpdateFields(obj)

        if obj.checkupdatedfields == 1

            obj.previousOptions = obj.options;

            % Help people based on typical vendor values. Set up so you
            % should still be able to enter custom value
            temp = obj.Prot.SpoilingCorrection.Mat;
            if strcmp(temp(1), 'Siemens')
                temp(1) = 50;
            elseif strcmp(temp(1), 'GE')
                temp(1) = 117;
            elseif strcmp(temp(1), 'Philips')
                temp(1) = 150;
            end
            
            obj.Prot.SpoilingCorrection.Mat = temp;
            

        elseif obj.checkupdatedfields == 2
         
            if obj.options.IncompleteSpoiling_RunSpoilingSimulations
                disp('Select directory to save simulations')
                obj.options.IncompleteSpoiling_OutputDirectory = uigetdir(pwd); 
                
                vfa_t1_spoil_hmri_corr_imperf_spoil(obj); 
    
            elseif obj.options.IncompleteSpoiling_LoadSpoilingSimulationsResults
 
                if ~obj.options.IncompleteSpoiling_RunSpoilingSimulations
                    disp('Load fitValues_SN.mat');
                    [FileName_SN, PathName_SN] = uigetfile('*.mat'); 
                    obj.fitValues_SN = load([PathName_SN filesep FileName_SN]); 
                else
                    disp('Select directory where you want values saved')
                    obj.options.IncompleteSpoiling_SeqSimDirectory = uigetdir(pwd);

                    disp('Load simulation fit values')
                    [FileName_sim,PathName_sim] = uigetfile('*.mat');
                         
                    obj.fitValues_sim = load([PathName_sim filesep FileName_sim]); 

                end
                
            end 

        end 
    end


    function FitResult = fit(obj,data)
         
        smallFlipApprox = obj.options.IncompleteSpoiling_smallAngleApprox;

        a1 = deg2rad(obj.Prot.VFAProtocol.Mat(1,1)); 
        a2 = deg2rad(obj.Prot.VFAProtocol.Mat(2,1));
        TR1 = deg2rad(obj.Prot.VFAProtocol.Mat(1,2)); 
        TR2 = deg2rad(obj.Prot.VFAProtocol.Mat(2,2)); 
        
        if ~isfield(data, 'Mask')
            data.Mask = zeros(size(data.LowFlipAngle));
            data.Mask(data.LowFlipAngle>100) =1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Apply spoiling correction:
        param.prot_name  = strcat('vfa_t2_', num2str(T2(3)));    % Used in output
        
        % I have wrote to a mat file
        simCoeffs = load(fullfile(param.outdir,[strrep(param.prot_name,' ',''),'_ABcoeff.mat']) );

        [T1corr, R1corr, S0corr] = vfa_t1_spol_calcMaps(data, a1, a2, TR1, TR2,...
                                                smallFlipApprox, simCoeffs);
          

        FitResult.T1corr = T1corr;
        FitResult.S0corr = S0corr;
        FitResult.R1corr = R1corr;
    end 
 
    end
end



























