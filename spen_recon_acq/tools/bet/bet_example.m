clear

data = load_nii_data('~/Dropbox/diffusion_toolbox/data/HARDI150.nii.gz');
bvals = load('~/Dropbox/diffusion_toolbox/data/HARDI150.bval');

% Get a brain mask
[data_brain, brain_mask] = bet('data',data, 'bvals',bvals);

data_slice = data(:,:,38,1);
brain_slice = data_brain(:,:,38,1);

multi = cat(3, rot90(data_slice), rot90(brain_slice));
figure(1)
montage(multi, 'DisplayRange', [0 2000]);
