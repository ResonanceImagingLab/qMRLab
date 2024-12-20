%% Code to make T1 and M0 maps for ihMT 

% Read the images 
[~, sp_mp2r_uni] = minc_read('spMP2RAGE_UNI_reg.mnc');
[~, sp_mp2r_inv2] = minc_read('spMP2RAGE_inv2_reg.mnc');
[~, b1] = minc_read('resampled_b1field.mnc'); 
[~, mask1] = minc_read('itkmask.mnc'); 

% Get b1 in proper image matrix orientation
b1 = permute(b1,[3 1 2]);

% View the images 
figure; imshow3Dfull(sp_mp2r_uni); 
figure; imshow3Dfull(sp_mp2r_inv2); 
figure; imshow3Dfull(b1); 
figure; imshow3Dfull(mask1);

% Flip b1 and mask 
b1 = flip(b1,1); 
b1 = flip(b1,2);
b1 = flip(b1, 3);
figure; imshow3Dfull(b1);

mask1 = flip(mask1, 1); 
mask1 = flip(mask1, 2);
mask1 = flip(mask1, 3); 
figure; imshow3Dfull(mask1);

% Write and save new b1 and mask 
[hdr, ~] = minc_read('spMP2RAGE_UNI_reg.mnc');
minc_write( 'b1.mnc', hdr, b1);
minc_write('mask.mnc', hdr, mask1); 

% Test it worked
[~, new_b1] = minc_read('b1.mnc'); figure; imshow3Dfull(new_b1); 
[~, new_mask] = minc_read('mask.mnc'); figure; imshow3Dfull(new_mask); 


%% T1 and M0 from the SPARSE MP2RAGE:

MP2RAGE.B0 = 3;           % in Tesla
MP2RAGE.TR = 5;           % MP2RAGE TR in seconds
MP2RAGE.TRFLASH = 6.4e-3; % TR of the GRE readout
MP2RAGE.TIs = [0.94 2.83];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices = [ceil(175/2) floor(175/2)];%  should be two values, number of excitations before k-space center, and number after. [Slices Per Slab * [PartialFourierInSlice-0.5  0.5] ]
MP2RAGE.FlipDegrees = [4 5];% Flip angle of the two readouts in degrees


MP2RAGEimg.img = sp_mp2r_uni; % load_untouch_nii(MP2RAGE.filenameUNI);
MP2RAGEINV2img.img = sp_mp2r_inv2; % load_untouch_nii(MP2RAGE.filenameINV2);
B1.img = b1;
brain.img = mask1;

tic
[ spT1map, spMP2RAGEcorrected, spAppmap2] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);
toc

spT1_map = spT1map.img;
spApp_mp2 = spAppmap2.img;

spT1_map = limitHandler(spT1_map, 0 , 6000);
spApp_mp2 = double(limitHandler(spApp_mp2, 0, 20000));

figure; imshow3Dfull(spT1_map, [300 2500],jet)
figure; imshow3Dfull(spApp_mp2 , [00 15000])

% Write T1 and M0 map 
minc_write('T1map.mnc', hdr, spT1_map); 
minc_write('M0map.mnc', hdr, spApp_mp2);

% Test they wrote properly 
[~, T1] = minc_read('T1map.mnc'); figure; imshow3Dfull(T1); 
[~, M0] = minc_read('M0map.mnc'); figure; imshow3Dfull(M0); 
