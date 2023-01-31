addpath(genpath('./'));
%%
CSparam=init;
CSparam.FT=1;
CSparam.Itnlim=20;
CSparam.TV=TVOP;
CSparam.XFM=1;
CSparam.xfmWeight=0;
% r2,4.3e1; r5,1.7e1; r1.5,1e2
CSparam.TVWeight= 4.3e1;      %5.3e1;

CSparam.data = image_recon_low;
% CSparam.data = im(:,:,13);
% CSparam.BW = 0.9*power(edge(:,:,13),10);
CSparam.BW = 0.8*power(edge_prior,1);
% CSparam.BW = ones(64,128);
% x0 = zeros(70,70); 
x0=zeros(res_sample*xy,res_sample);

for n=1:1:8
    x1=CG1(x0,CSparam);
    x0=x1;
end

figure();imshowMRI(abs(image_recon_low),[0,0.6*max(abs(image_recon_low(:)))]);
figure();imshowMRI(abs(x1),[0,0.6*max(abs(x1(:)))]);
%% zoom in 1
figure(101);imshowMRI(abs(imresize(image_recon_low(56:90,61:90),[105,90])),...
    [min(abs(image_recon_low),[],'all'),0.5*max(abs(image_recon_low),[],'all')]);
figure(102);imshowMRI(abs(imresize(x1(56:90,61:90),[105,90])),...
    [min(abs(x1),[],'all'),0.5*max(abs(x1),[],'all')]);
%% zoom in 2
figure(101);imshowMRI(abs(imresize(image_recon_low(31:65,31:60),[105,90])),...
    [min(abs(image_recon_low),[],'all'),0.5*max(abs(image_recon_low),[],'all')]);
figure(102);imshowMRI(abs(imresize(x1(31:65,31:60),[105,90])),...
    [min(abs(x1),[],'all'),0.5*max(abs(x1),[],'all')]);
