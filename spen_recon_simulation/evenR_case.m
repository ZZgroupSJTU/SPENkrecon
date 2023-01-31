%% declaration
% it's a k space convertion 
% large bandwidth subsampled SPEN to narrow bandwidth full-sampled EPI
% valid resolution remian
%% setting
addpath(genpath('./'));
close all;clc;clear;
RES = 8192;
FOV = 30;
res_sample = 128;
R = 2;
Q = res_sample *R /2 /FOV^2;
res_full = res_sample *R;
sampling_pattern = 1:R:res_full;

%% generate data
% Phantom = phantom('Modified Shepp-Logan',RES);
% Phantom_noised = imnoise(Phantom,'gaussian'); xy = 1;
load('phantom_brain_256x320x320.mat'); xy=1.25;
im0 = rot90(phantom_brain_256x320x320(:,:,180));
% im0 = phantom_brain_256x320x320(:,:,230);
figure(99);imshowMRI(im0,[]);
Phantom_noised = imresize(im0,[RES*xy,RES]);

full_discrete_FOV = (-RES/2:RES/2-1)/RES*FOV;
quadratic = exp(1i*2*pi*Q*full_discrete_FOV.^2);
Phantom_noised_quadratic_phase = bsxfun(@times,Phantom_noised,quadratic);
ksp_EPI_full = cfftn(Phantom_noised, [1,2]);
ksp_SPEN_full = cfftn(Phantom_noised_quadratic_phase, [1,2]);

% ksp_EPI = crop(ksp_EPI_full,[res_full*xy/R, res_full]);
ksp_EPI = crop(ksp_EPI_full,[res_full*xy/R, res_full]);
ksp_SPEN = crop(ksp_SPEN_full,[res_full*xy/R, res_full]);
ksp_SPEN_sample = ksp_SPEN(:,sampling_pattern);

%% low frequency recon
if mod(R,2)~=1
    res_pad = res_sample *(R+1);
    ksp_EPI = crop(ksp_EPI_full,[res_full*xy/R, res_pad]);
    %ksp_EPI = crop(ksp_EPI_full,[res_pad*xy, res_pad]);
    R_process = R+1;
else
    res_pad = res_full;
    R_process = R;
end

k_in = (-res_pad/2:res_pad/2-1)/FOV;
k_out = (-res_full/2:res_full/2-1)/FOV;
A = A_k2k(Q,k_in,k_out);
Asub = A(sampling_pattern,:);
Asublow = crop(Asub,[res_sample,res_sample]);
k_recon_low = (Asublow\ksp_SPEN_sample.').';
image_recon_low = cifftn(k_recon_low,[1,2])/R;

figure(1);imshowMRI(abs(image_recon_low));

% figure(101);imshowMRI(abs(imresize(image_recon_low(96:130,51:80),[105,90])),...
%     [min(abs(image_recon_low),[],'all'),max(abs(image_recon_low),[],'all')]);
%% extract edge ghost from k-space high frequency regions
% Asub_parts = zeros(res_sample,res_sample,R);
% kepi_parts = zeros(res_sample,res_sample,R);
% Ps = zeros(res_sample,res_sample,R);
for r = 1:R_process
    if mod(R_process,2)==1
        Asub_parts(:,:,r) = Asub(:,res_sample*(r-1)+1:res_sample*r);
        A_parts(:,:,r) = A(:,res_sample*(r-1)+1:res_sample*r);
        kepi_parts(:,:,r) = ksp_EPI(:,res_sample*(r-1)+1:res_sample*r);
    else
        if r < R_process
            Asub_parts(:,:,r) = Asub(:,res_sample*(r-1/2)+1:res_sample*(r+1/2));
            kepi_parts(:,:,r) = ksp_EPI(:,res_sample*(r-1/2)+1:res_sample*(r+1/2));
        else
            Asub_parts(:,:,r) = cat(2,Asub(:,1:res_sample/2),Asub(:,res_sample*(R_process-1/2)+1:end));
            kepi_parts(:,:,r) = cat(2,ksp_EPI(:,1:res_sample/2),ksp_EPI(:,res_sample*(R_process-1/2)+1:end));
        end
    end
    Ps(:,:,r) = Asublow\Asub_parts(:,:,r);
    quadratic_k(:,:,r) = (Asub_parts(:,:,r)*kepi_parts(:,:,r).').';
    quadratic_k_all(:,:,r) = (A_parts(:,:,r)*kepi_parts(:,:,r).').';
    Psmodkepi_parts(:,:,r) = (Ps(:,:,r)*kepi_parts(:,:,r).').';
end
Psmodim = cifftn(Psmodkepi_parts,[1,2]);
quadratic_im = cifftn(quadratic_k,[1,2]);
% m = ones(1,1,R_process); m(1,1,:) = [5,1,5];
m = 1;
im = cifftn(kepi_parts,[1,2]);
figure(2);imshowMRI(abs(quadratic_im.*m),[0,0.6*max(abs(quadratic_im(:)))],[R_process,1]);
figure(3);imshowMRI(abs(Psmodim.*m),[0,0.6*max(abs(Psmodim(:)))],[R_process,1]);
figure(4);imshowMRI(abs(im.*m),[0,0.6*max(abs(im(:)))],[R_process,1]);
%% extract edge ghost detection (differential)
diff = eye(res_sample); diff = eye(res_sample)-[diff(:,2:end),diff(:,1)];
for r = 1:R_process
    if floor((R_process+1)/2)==r
        continue;
    else
        ktemp = (Ps(:,:,r)*k_recon_low.').';
        k_ps_low(:,:,r) = ktemp;
        imtemp = cifftn(ktemp,[1,2]);
        edge_parts(:,:,r) = (diff*imtemp.').';
    end
end
edge = mean(abs(edge_parts),3);
edge_prior = edge/max(edge,[],'all');
figure(5);imshowMRI(abs(edge_prior),[]);

%% extract edge ghost detection (lowpass filter)
h = fspecial('gaussian',res_sample*xy,13);
h = crop(h,[res_sample*xy,res_sample]);
h = h/max(h,[],'all')*0.95;
hx = 1-h;
for r = 1:R_process
    if floor((R_process+1)/2)==r
        continue;
    else
        ktemp = (Ps(:,:,r)*k_recon_low.').';
        k_ps_low(:,:,r) = ktemp;
        ktemp = ktemp .* hx;
        edge_parts(:,:,r) = cifftn(ktemp,[1,2]);
    end
end
edge = mean(abs(edge_parts),3);
edge_prior2 = edge/max(edge,[],'all');
figure(6);imshowMRI(abs(edge_prior2),[]);