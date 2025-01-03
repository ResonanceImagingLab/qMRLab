function sig = ihMT_generate_BSF_scaling_v1( inputMag, Params, outputSamplingTable , gm_m, fft_gm_m)

% Over long excitation trains, you have a change in the magnetization.
% if you have a map that is being calculated off of the first echo in the
% train, it is assuming the signal == across the train.

% gm_m is the grey matter segmentation that is the same reconstructred size as the sampling table
% fft_gm_m is the k-space values of the above (complex!). It is expected
% that fftshift has already been run on this value (middle is the
% brightest)


%% Fill table with echo train values
outKspace_s = ihMT_fillKspaceSamplingTable_v2( inputMag, outputSamplingTable, Params);

%% Interpolate missing grappa lines
outKspace_s = ihMT_interpolateMissingGrappaLines( outKspace_s);

%% Assume constant values in readout direction 
sim3d_m = repmat(outKspace_s, [1,1,Params.Slices]);
targetSize = size(fft_gm_m);
numSlices = size(sim3d_m, 3);
sim3d_m_resized = zeros(targetSize);
    for slice = 1:min(numSlices, targetSize(3))
        sim3d_m_resized(:,:,slice) = imresize(sim3d_m(:,:,slice), [targetSize(1), targetSize(2)]);
    end

%% Scale the k-space by brain weighting to get 'Brain-spread-function'
bsf = sim3d_m_resized .* fft_gm_m;
%bsf = sim3d_m .* fft_gm_m;

b_v = abs(ifftn(ifftshift(bsf)));

sig = mean( b_v( gm_m>0 ));



% Code for viewing/debugging
% figure; imagesc( outKspace_s);
% axis image; colorbar;