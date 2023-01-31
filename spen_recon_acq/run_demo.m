addpath(genpath('./'));
close all;clc;clear;
%% load and set
% load MAT,para,
% set prior power,weight

% load('./Data/ACR_3x3_Q105_R3.mat');
% pow = 12;
% weight = 9e1;

load('./Data/Brain_Q105.mat');
pow = 10;
weight = 5.3e1;

% load('./Data/Orange_1p5x1p5_Q96_R3.mat');
% pow = 13;
% weight = 8e1;

[imw,im,edgew,edge] = SPEN_recon_SJTU_mat(mat,SPENpara);
%% optimization
close all;
edge0 = 0.9*power(edge,pow); 
imopt = opt_process(im,edge0,weight);
figure();imshowMRI(abs(edge));
%% compare
imshape = size(imopt);
Npe = imshape(1);
Nro = imshape(2);
winf=(hann(Npe)*hann(Nro).').^(1/3);
imopt_w=imresize(cifftn(bsxfun(@times,cfftn(imopt,[1 2]),winf),[1 2]),[Npe*2,Nro*2],'lanczos2');
im0_w=imresize(cifftn(bsxfun(@times,cfftn(im,[1 2]),winf),[1 2]),[Npe*2,Nro*2],'lanczos2');
edge_w=imresize(cifftn(bsxfun(@times,cfftn(edge0,[1 2]),winf),[1 2]),[Npe*2,Nro*2],'lanczos2');
ma = max(abs(imopt_w),[],'all');
figure();imshow(abs(imopt_w-im0_w),0.08*ma*[0,0.4]);title('error between Approximate & Optimized');
figure();imshow(power(abs(imopt_w),0.65),[]);title('Optimized recon, rescaled');
figure();imshow(power(abs(im0_w),0.65),[]);title('Approximate recon, rescaled');
%% ACR local
% figure();imshow(abs(imresize(imopt_w(56:90,36:70),[70,70])),[ma*0.2,ma*0.9]);title('Zoomed in, Optimized recon, rescaled');
% figure();imshow(abs(imresize(im0_w(56:90,36:70),[70,70])),[ma*0.2,ma*0.9]);title('Zoomed in, Approximate recon, rescaled');
%% Brain local
m1 = power(abs(imopt_w),0.65);
m0 = power(abs(im0_w),0.65);
M = max([m0,m1],[],'all');
% Q56
% figure();imshow(abs(imresize(m1(20:54,64:98),[70,70])),M*[0.08,0.5]);title('Zoomed in, Optimized recon, rescaled');
% figure();imshow(abs(imresize(m0(20:54,64:98),[70,70])),M*[0.08,0.5]);title('Zoomed in, Approximate recon, rescaled');
% Q70
% figure();imshow(abs(imresize(m1(33:57,30:64),[70,70])),M*[0.08,0.5]);title('Zoomed in, Optimized recon, rescaled');
% figure();imshow(abs(imresize(m0(33:57,30:64),[70,70])),M*[0.08,0.5]);title('Zoomed in, Approximate recon, rescaled');
% Q105
figure();imshow(abs(imresize(m1(51:85,71:105),[70,70])),M*[0.08,0.5]);title('Zoomed in, Optimized recon, rescaled');
figure();imshow(abs(imresize(m0(51:85,71:105),[70,70])),M*[0.08,0.5]);title('Zoomed in, Approximate recon, rescaled');
