function [ T1corr, MP2RAGEcorrected, M0map_corr] = ihMT_T1B1correctpackageTFL_withM0( B1img, MP2RAGEimg, MP2RAGEINV2nii, MP2RAGE, brain, InvEff)
% Adapted from T1B1correctpackageTFL, to provide a B1 corrected value for M0
% as well.

% While this outputs a B1 corrected UNI image, you need two images for the
% M0 image. So the correction needs to be done further up. 

% usage
%
% [  T1corr, MP2RAGEcorrected, M0map_corr] = CR_T1B1correctpackageTFL_withM0(B1, MP2RAGEimg,T1,MP2RAGE,brain,varargin)
%

% B1 and MP2RAGEimg (and T1img) are the nii structures resulting from loading
% the MP2RAGE (MP2RAGEimg, T1img) and the result of some B1 mapping technique
% with load_nii or load_untouch_nii
%
% The variable B1 is compulsory
%
% Only MP2RAGEimg or the T1img have to be loaded (I usually use the MP2RAGEimg)
%
% The variables that are not loaded can be simply left empty []
%
% MP2RAGE variable contains all the relevant sequence
% information as delailed below
%
% B1 will be given in relative units => 1 if it was correct; values can vary from 0-2
%
%     MP2RAGE.B0          = 7;                  % In Tesla
%     MP2RAGE.TR          = 6;                  % MP2RAGE TR in seconds
%     MP2RAGE.TRFLASH     = 6.7e-3;             % TR of the GRE readout
%     MP2RAGE.TIs         = [800e-3 2700e-3];   % Inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
%     MP2RAGE.NZslices    = [40 80];            % Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
%     MP2RAGE.FlipDegrees = [4 5];              % Flip angle of the two readouts in degrees
%
% Brain can be an image in the same space as the MP2RAGE that has zeros
% where there is no need to do any T1/B1 calculation (can be a binary mask
% or not). if left empty the calculation is done everywhere
%
% Additionally the inversion efficiency of the adiabatic inversion can be
% set as a last optional variable. Ideally it should be 1.
% In the first implementation of the MP2RAGE the inversino efficiency was
% measured to be ~0.96
%
% Outputs are:
%  T1corr  - T1map corrected for B1 bias
%  MP2RAGEcorr - MP2RAGE UNI image corrected for B1 bias
%
% Please cite:
%  Marques, J.P., Gruetter, R., 2013. New Developments and Applications of the MP2RAGE Sequence - Focusing the Contrast and High Spatial Resolution R1 Mapping. PLoS ONE 8. doi:10.1371/journal.pone.0069294
%  Marques, J.P., Kober, T., Krueger, G., van der Zwaag, W., Van de Moortele, P.-F., Gruetter, R., 2010a. MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. NeuroImage 49, 1271�1281. doi:10.1016/j.neuroimage.2009.10.002

% 
% B1img = B1;
% MP2RAGEimg = MP2RAGEimg;
% MP2RAGEINV2nii = MP2RAGEINV2img;
% MP2RAGE = MP2RAGE;
% brain = brain;
% InvEff = 0.96;


%% Parse the input arguments

if nargin < 6 || isempty(InvEff)
    InvEff = 0.99;
end

if nargin < 5 || isempty(brain)
    if isempty(MP2RAGEimg)
        brain = T1img;
    else
        brain = MP2RAGEimg;
    end
    brain.img = ones(size(brain.img));
end


%% definition of range of B1s and T1s and creation of MP2RAGE lookupvector to make sure the input data for the rest of the code is the MP2RAGEimg

% [MP2RAGE.Intensity, MP2RAGE.T1vector] = MP2RAGE_lookuptable(2, MP2RAGE.TR, MP2RAGE.TIs, MP2RAGE.FlipDegrees, MP2RAGE.NZslices, MP2RAGE.TRFLASH, 'normal', InvEff);

% C.R. removed and replaced the T1 input with MP2RAGEINV2nii. 
% if isempty(MP2RAGEimg)
%     T1img.img      = double(T1img.img)/1000; % C.R. notes that it is assuming here that your T1 map is in milliseconds
%     MP2RAGEimg.img = reshape(interp1(MP2RAGE.T1vector, MP2RAGE.Intensity, T1img.img(:)), size(B1img.img));    
%     % MP2RAGEimg.img(isnan(MP2RAGEimg.img)) = -0.5; C.R. was this...
%     MP2RAGEimg.img(isnan(MP2RAGEimg.img)) = 0.5;
% else
%     MP2RAGEimg.img = double(MP2RAGEimg.img)/4095-0.5;
% end

MP2RAGEimg.img = double(MP2RAGEimg.img)/4095-0.5;
%% now the fun starts

% creates a lookup table of MP2RAGE intensities as a function of B1 and T1

B1_vector = 0.005:0.05:1.9;
T1_vector = 0.05:0.05:5; % C.R. changed to match the one in MP2RAGE_lookuptable() -> could just crop out the values from MP2RAGE LOOKUPTABLE

MP2RAGEmatrix = zeros(length(B1_vector), length(T1_vector));
Inv2Sig = zeros(length(B1_vector), length(T1_vector));

for k = 1:length(B1_vector)
    
    b1val = B1_vector(k);
    
    % Calculate the lookup table for this B1 value
    [Intensity, T1vector, uncombIntensity] = MP2RAGE_lookuptable(2, MP2RAGE.TR, MP2RAGE.TIs, b1val*MP2RAGE.FlipDegrees, MP2RAGE.NZslices, MP2RAGE.TRFLASH, 'normal', InvEff,1);
    
    % Determine the MP2RAGE intensity for a given T1 and B1map value 
    MP2RAGEmatrix(k,:) = interp1(T1vector, Intensity, T1_vector);
    
    % Store the 2nd Inversion signal values for later:
    Inv2Sig(k,:) = uncombIntensity(:,2)';
end


% C.R. view the surface result
% [b1mesh, t1mesh] = meshgrid(B1_vector, T1_vector);
% figure; surf( b1mesh', t1mesh', MP2RAGEmatrix);


%% make the matrix MP2RAGEMatrix into T1_matrix(B1, ratio)

% C.R. increasing the length of this vector removed NaN's
MP2RAGE_vector = linspace(-0.5, 0.5, 150);


for k = 1:length(B1_vector)
    %b1val = B1_vector(k);
    try
        T1matrix(k,:) = interp1(MP2RAGEmatrix(k,:), T1_vector, MP2RAGE_vector, 'pchirp'); 
    catch
        temp              = MP2RAGEmatrix(k,:); 
        temp(isnan(temp)) = linspace(-0.5-eps, -1, sum(isnan(temp(:))));
        temp              = interp1(temp, T1_vector, MP2RAGE_vector);
        
        T1matrix(k,:) = temp;       
    end
    
end


%% correcting the estimates of T1 and B1 iteratively

T1temp                        = MP2RAGEimg;

brain.img(B1img.img==0)       = 0;
brain.img(MP2RAGEimg.img==0)  = 0;
T1temp.img(brain.img==0)      = 0;
T1temp.img(brain.img==1)      = 0;
B1img.img(brain.img==0)       = 0;

% temp1                         = squeeze(T1temp.img(:, end/2, :));

T1temp.img(brain.img~=0)      = interp2(MP2RAGE_vector, B1_vector, T1matrix, MP2RAGEimg.img(brain.img~=0), B1img.img(brain.img~=0));
T1temp.img(isnan(T1temp.img)) = 4;  % Set NaN to 4sec: When T1s are very long, you can get some nan out of the lookup table (that happens for some protocols for CSF)

% temp2                         = squeeze(T1temp.img(:, end/2, :));




%% creates an MP2RAGEcorrected image and puts both the B1 and T1 in the ms scale
if nargout > 1
    [MP2RAGE.Intensity, MP2RAGE.T1vector] = MP2RAGE_lookuptable(2, MP2RAGE.TR, MP2RAGE.TIs, MP2RAGE.FlipDegrees, MP2RAGE.NZslices, MP2RAGE.TRFLASH, 'normal', InvEff);
    
    MP2RAGEcorrected     = MP2RAGEimg;
    MP2RAGEcorrected.img = reshape(interp1(MP2RAGE.T1vector, MP2RAGE.Intensity, T1temp.img(:)), size(T1temp.img));
    MP2RAGEcorrected.img(isnan(MP2RAGEcorrected.img)) =- 0.5;
    MP2RAGEcorrected.img = round(4095*(MP2RAGEcorrected.img + 0.5));
end

T1temp.img = (T1temp.img)*1000;         % Retain precision when saved as integer
T1corr = T1temp;


%% At this point, we have the T1 values calculated and corrected for B1.
% Need to determine the M0 value, provided the value of INV2 image,
% calculated T1, and acquired B1 value. 
M0map_corr = T1corr; % copy the header
M0map_corr.img = zeros(size(T1corr.img)); % clear the image 

%% Convert to vectors for easier fitting

% find indices of valid voxels
q = find( (brain.img(:) > 0));

t1_l = T1temp.img(q)./1000; % convert back to seconds...
b1_l = B1img.img(q);

%% Solve the signal matrix as a function of T1 and B1

[x, y] = ndgrid(B1_vector, T1_vector);
F = griddedInterpolant(x,y,Inv2Sig);

% Use the interpolant to make an image of relative signal provided the B1
% and T1 value
m0_l = F(b1_l,t1_l);

% The M0 map will be the 2nd inversion image divided by this map
M0map_corr.img(q) = m0_l;     % fill the image 
M0map_corr.img = MP2RAGEINV2nii.img ./ M0map_corr.img  .* brain.img;     % fill the image


% figure; imshow3Dfull(MP2RAGEINV2nii.img ./M0map_corr.img .* brain.img )

% for i = 1:length(t1_l)
%     
%     % Extract the values based on b1 measured
%     vect1 = interp1 (B1_vector , Inv2Sig_mat, b1_l(i));
%     
%     % next extract values based on T1 measured
%     vect2 = interp1 (T1_vector , squeeze(vect1), t1_l(i));
%     
%     % Now use the measured signal to extract an M0 value
%     x = interp1 (vect2, M0_vector, img_l(i))
% 
% end





%% C.R. I moved all the figure stuff to the bottom.
showimages = 0;
if showimages==1    
    imagesc(MP2RAGE_vector, B1_vector, T1matrix, [0.4 5])
    colorbar
    xlabel('MP2RAGE', 'FontSize',12, 'FontWeight','bold')
    ylabel('B_1', 'FontSize',12, 'FontWeight','bold')
    title('T_1 look-up table', 'FontSize',12, 'FontWeight','bold')
    set(gca, 'FontSize',12, 'LineWidth',2)
end



%% sanity check to see how B1 sensitive your sequence was

% H1 = figure(1);
% set(H1, 'Color',[1 1 1], 'Name','B1-sensitivity');
% hold off
% 
% for B1val = 0.6:0.2:1.4
%     
%     [MP2RAGEamp, T1vector] = MP2RAGE_lookuptable(2, MP2RAGE.TR, MP2RAGE.TIs, B1val*MP2RAGE.FlipDegrees, MP2RAGE.NZslices, MP2RAGE.TRFLASH, 'normal');
%     
%     plot(MP2RAGEamp, T1vector, 'color',[0.5 0.5 0.5]*B1val, 'Linewidth',2)
%     hold on
%     
% end
% legend('B1=0.6', 'B1=0.8', 'B1=1', 'B1=1.2', 'B1=1.4')
% 
% % examples of T1 values at 3T
% 
% if ~isfield(MP2RAGE,'B0')
%     T1WM  = 1.1;
%     T1GM  = 1.85;
%     T1CSF = 3.5;
% elseif MP2RAGE.B0==3
%     T1WM  = 0.85;
%     T1GM  = 1.35;
%     T1CSF = 2.8;
% else
%     % examples of T1 values at 7T
%     T1WM  = 1.1;
%     T1GM  = 1.85;
%     T1CSF = 3.5;
% end
% 
% plot([-0.5 0.5], [T1CSF T1CSF; T1GM T1GM; T1WM T1WM]', 'Linewidth',2)
% text(0.35, T1WM,  'White Matter')
% text(0.35, T1GM,  'Grey Matter')
% text(0.35, T1CSF, 'CSF')
% ylabel('T1');
% xlabel('MP2RAGE');

% H2 = figure(2);
% set(H2, 'Color',[1 1 1], 'Name','T1-correction');
% imagesc(temp2 - temp1)
% colorbar
% title('T1 correction');
% colormap(gray)