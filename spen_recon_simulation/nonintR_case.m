%% declaration
% it's a k space convertion 
% large bandwidth subsampled SPEN to narrow bandwidth full-sampled EPI
% valid resolution remian
%% setting
addpath(genpath('./'));
close all;clc;clear;
RES = 4096;
FOV = 30;
res_sample = 128;
Numerator = 3; Denumerator = 2; R = Numerator/Denumerator;
Q = res_sample *R /2 /FOV^2;
res_full = res_sample *R;
sampling_pattern = 1:Numerator:res_full*Denumerator;

%% generate data
% Phantom = phantom('Modified Shepp-Logan',RES);
% Phantom_noised = imnoise(Phantom,'gaussian'); xy = 1;
load('phantom_brain_256x320x320.mat'); xy=1.25;
im0 = rot90(phantom_brain_256x320x320(:,:,180));
figure(99);imshowMRI(im0,[]);

% EPI
Phantom_noised = imresize(im0,[RES*xy,RES]);
ksp_EPI_full = cfftn(Phantom_noised, [1,2]);
ksp_EPI = crop(ksp_EPI_full,[res_full*xy/R, res_full]);
% SPEN R
Phantom_noised_R = cat(2,zeros(RES*xy,RES*((Denumerator-1)/2)),...
    Phantom_noised,zeros(RES*xy,RES*((Denumerator-1)/2)));
RES_R = RES*Denumerator; FOV_R = FOV*Denumerator;
full_discrete_FOV_R = (-RES_R/2:RES_R/2-1)/RES_R*FOV_R;
quadratic = exp(1i*2*pi*Q*full_discrete_FOV_R.^2);
Phantom_noised_quadratic_phase_R = bsxfun(@times,Phantom_noised_R,quadratic);

ksp_SPEN_full = cfftn(Phantom_noised_quadratic_phase_R, [1,2]);
ksp_SPEN = crop(ksp_SPEN_full,[res_full*xy/R, res_full*Denumerator]);
ksp_SPEN_sample = ksp_SPEN(:,sampling_pattern);

%% low frequency recon
if mod(R,2)~=1
    R_process = 2*floor((R+1)/2)+1;
    res_pad = res_sample * R_process;
    ksp_EPI = crop(ksp_EPI_full,[res_full*xy/R, res_pad]);
else
    res_pad = res_full;
    R_process = R;
end

k_in = (-res_pad/2:res_pad/2-1)/FOV;
k_out = (-res_full/2:R:res_full/2-1)/FOV;
Asub = A_k2k(Q,k_in,k_out);
Asublow = crop(Asub,[res_sample,res_sample]);
k_recon_low = (Asublow\ksp_SPEN_sample.').';
image_recon_low = cifftn(k_recon_low,[1,2]);

figure(1);imshowMRI(abs(image_recon_low),[]);
%% extract edge ghost from k-space high frequency regions
% Asub_parts = zeros(res_sample,res_sample,R);
% kepi_parts = zeros(res_sample,res_sample,R);
% Ps = zeros(res_sample,res_sample,R);
for r = 1:R_process
    if mod(R_process,2)==1
        Asub_parts(:,:,r) = Asub(:,res_sample*(r-1)+1:res_sample*r);
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
    Psmodkepi_parts(:,:,r) = (Ps(:,:,r)*kepi_parts(:,:,r).').';
end
Psmodim = cifftn(Psmodkepi_parts,[1,2]);
quadratic_im = cifftn(quadratic_k,[1,2]);
m = ones(1,1,R_process); m(1,1,:) = [5,1,5];
im = cifftn(kepi_parts,[1,2]);
figure(2);imshowMRI(abs(quadratic_im.*m),[],[R_process,1]);
figure(3);imshowMRI(abs(Psmodim.*m),[],[R_process,1]);
figure(4);imshowMRI(abs(im.*m),[],[R_process,1]);
%% extract edge ghost detection
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
