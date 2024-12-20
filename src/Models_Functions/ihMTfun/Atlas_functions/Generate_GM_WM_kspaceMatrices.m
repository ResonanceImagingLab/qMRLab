%% Generate_GM_WM_kspaceMatrices
% This is mainly needed for centric-out encoding on Siemens.

DATADIR = '/Users/amiedemmans/Documents/GitHub/OptimizeIHMTimaging/kspaceWeighting/Atlas_reference/mni_icbm152_nlin_sym_09a_minc2/';
OutputDir =  '/Users/amiedemmans/Documents/ihMT_images_from_Chris/atlas_images/';



% Load the maps
mtw_fn = {'mni_icbm152_csf_tal_nlin_sym_09a.mnc';'mni_icbm152_gm_tal_nlin_sym_09a.mnc'; 'mni_icbm152_wm_tal_nlin_sym_09a.mnc';'DeepStructureMask.mnc'};
                                                
% load images
%comb_mtw = zeros(224, 256, 176, 18);

for i = 1:length(mtw_fn)
    fn = strcat(DATADIR,mtw_fn{i});
    [hdr, img] = minc_read(fn);
    comb_mtw(:,:,:,i) = img; %.img;
end

[temp, temp2] = niak_read_vol('/media/chris/data8tb/Research/MagTransfer/20240806_pediatricIHMT/minc/ped_ihmt_20240806_085237_15d1_mri.mnc.gz');



figure; imshow3Dfull(comb_mtw(:,:,:,1)  )
figure; imshow3Dfull( temp2 )


Params.Orientation = 'Axial';

switch Params.Orientation
    case  'Axial'
        comb_mtw2 = comb_mtw; % already in this 
    case 'Sagittal'
        comb_mtw2 = permute(comb_mtw, [2,3,1,4]);
    case 'Coronal'
        comb_mtw2 = permute(comb_mtw, [3,1,2,4]);
    otherwise 
        disp( 'Set Params.Orientation to either Axial, Sagittal, or Coronal')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.7 mm
%wanted_size=[128,140,88,4];
wanted_size = [192,256,256,4];

reduc = wanted_size./size(comb_mtw2);
% to avoid scaling issues, first 3 values should be the same, remove rows
% and columns from comb_mtw
minScale = max(reduc(1:3));

% find matrix size differences
[x,y,z] = size(comb_mtw2(:,:,:,1));
matDif = round([x*minScale - wanted_size(1), ...
                y*minScale - wanted_size(2), ...
                z*minScale - wanted_size(3)]);

matDif = matDif + mod(matDif,2); % force even

% Remove extra
comb_resz = comb_mtw2;
if matDif(1) > 0
    comb_resz(1:matDif(1)/2, :, :, :) = [];
    comb_resz(end-matDif(1)/2+1:end,:,:, :) = [];
end

if matDif(2) > 0
    comb_resz( :, 1:matDif(2)/2, :, :) = [];
    comb_resz( :, end-matDif(2)/2+1:end,:, :) = [];

end

if matDif(3) > 0
    comb_resz( :, :, 1:matDif(3)/2, :) = [];
    comb_resz(:, :, end-matDif(3)/2+1:end, :) = [];
end

reduc = (wanted_size./size(comb_resz));

%comb_mtw3=imresizen( comb_resz, wanted_size./size(comb_resz));
comb_mtw3=atlas_imresizen( comb_resz, reduc);

[x,y,z] = size(comb_mtw3(:,:,:,1));

figure; imshow3Dfull(comb_mtw3(:,:,:,1)  )

%% Brain tissue:
brain_m = zeros(x,y,z);
brain_m(comb_mtw3(:,:,:,2) >=0.05) = 1; % add extra, then remove
brain_m(comb_mtw3(:,:,:,3) >=0.05) = 1;
fft_brain_m = fftshift(fftn(brain_m));
figure; imshow3Dfull( brain_m )

%% Save results.

save( strcat(OutputDir,'GM_seg_MNI_152_image.mat'), 'brain_m')
save( strcat(OutputDir,'GM_seg_MNI_152_kspace.mat'), 'fft_brain_m')


%% 1.2 mm

wanted_size=[140,182,128,4];

reduc = wanted_size./size(comb_mtw2);
% to avoid scaling issues, first 3 values should be the same, remove rows
% and columns from comb_mtw
minScale = max(reduc(1:3));

% find matrix size differences
[x,y,z] = size(comb_mtw2(:,:,:,1));
matDif = round([x*minScale - wanted_size(1), ...
                y*minScale - wanted_size(2), ...
                z*minScale - wanted_size(3)]);

matDif = matDif + mod(matDif,2); % force even

% Remove extra
comb_resz = comb_mtw2;
if matDif(1) > 0
    comb_resz(1:matDif(1)/2, :, :, :) = [];
    comb_resz(end-matDif(1)/2+1:end,:,:, :) = [];
end

if matDif(2) > 0
    comb_resz( :, 1:matDif(2)/2, :, :) = [];
    comb_resz( :, end-matDif(2)/2+1:end,:, :) = [];

end

if matDif(3) > 0
    comb_resz( :, :, 1:matDif(3)/2, :) = [];
    comb_resz(:, :, end-matDif(3)/2+1:end, :) = [];
end

reduc = wanted_size./size(comb_resz);


comb_mtw3=atlas_imresizen( comb_resz, wanted_size./size(comb_resz));

[x,y,z] = size(comb_mtw3(:,:,:,1));

figure; imshow3Dfull(comb_mtw3(:,:,:,1)  )

%% Brain tissue:
brain_m = zeros(x,y,z);
brain_m(comb_mtw3(:,:,:,2) >=0.05) = 1; % add extra, then remove
brain_m(comb_mtw3(:,:,:,3) >=0.05) = 1;
fft_brain_m = fftshift(fftn(brain_m));
figure; imshow3Dfull( brain_m )
%% Save results.

save( strcat(OutputDir,'GM_seg_MNI_152_image1p2.mat'), 'brain_m')
save( strcat(OutputDir,'GM_seg_MNI_152_kspace1p2.mat'), 'fft_brain_m')



