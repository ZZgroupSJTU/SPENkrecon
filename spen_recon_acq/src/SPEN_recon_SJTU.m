
function [im_final_w,im_final]=SPEN_recon_SJTU(datfile,add1)
% Simple SPEN reconstruction pipeline for multi-shot SPEN acquisition (only even-shot acquisiton without parallel acceleration is tested)
% The reconstruction corrects shot-shot phase variations by comparing the phases of shot-based SPEN reconstructed images
% No iteration is needed for this pipeline
% (c) Zhiyong Zhang@SJTU 2021.10.24
%% Prepare the dataset
disp('reading the rawdata from siemens twix file');
twix_obj_raw=mapVBVD(datfile,'removeOS','ignoreSeg','rampSampRegrid');
ksp_raw=SPEN_parse_siemenstwix(twix_obj_raw);
[Nx,Nc,Netl,Ns,Nseg,Nb]=size(ksp_raw);
ksp_spen=permute(reshape(permute(ksp_raw,[1 5 3 2 4 6]),[Nx Netl*Nseg Nc Ns Nb]),[2 1 3 4 5 6]);
disp('reading the rawdata from siemens twix file, done');
%% Get SPEN parameter
disp('Extracting SPEN reconstruction parameters');
SPENpara=getSPENreconParaFromTwixobj(twix_obj_raw{end},'SPEN_diff');
Nacq=size(ksp_spen,1);
Nfull=2*SPENpara.recon.Q;
Q=-SPENpara.recon.Q/(SPENpara.infor.PEfov)^2;
% k_source=(-Nfull/2:Nfull/2-1)*1/SPENpara.infor.PEfov + 0*SPENpara.recon.PEshiftnorm;
k_source=(-Nacq*3/2:Nacq*3/2-1)*1/SPENpara.infor.PEfov + 0*SPENpara.recon.PEshiftnorm;
k_response=(-Nacq/2:Nacq/2-1)*Nfull/Nacq*1/SPENpara.infor.PEfov ;
A= A_k2k(Q,k_source,k_response);
disp('Extracting SPEN reconstruction parameters, done');
%% SPEN reconstruction
disp('Doing SPEN reconstruction');
Npe=size(ksp_spen,1);
A_ref=crop(A(1:Nseg:end,:),[Npe/Nseg Npe/Nseg]);
AHA=A_ref'*A_ref;
[~,S,~]=svd(AHA);
lambda=0.4*S(1,1);
AH_ref=pinv(AHA + lambda *eye(size(AHA)))*A_ref';   
dsize=size(ksp_spen);
dsize(1)=dsize(1)/Nseg;
ksp_ref=reshape(AH_ref*ksp_spen(1:Nseg:end,:),dsize);
im_ref=cifftn(ksp_ref,[1 2]);

ksp_spen_c=ksp_spen;
for m=2:Nseg
    A_currentseg=crop(A(m:Nseg:end,:),[Npe/Nseg Npe/Nseg]);
    ksp_currentseg=reshape(pinv(A_currentseg)*ksp_spen(m:Nseg:end,:),dsize);
    im_currentseg=cifftn(ksp_currentseg,[1 2]);
    im_phc=bsxfun(@times,sos(im_currentseg,3),exp(1i*angle(im_currentseg.*conj(im_ref))));
    ph_c=-angle(cifftn(bsxfun(@times,cfftn(im_phc,[1 2]),hamming(Npe/Nseg)*(hamming(Nx).').^(2)),[1 2]));
%     ph_c=repmat(reshape(ph_c,[Npe/Nseg,Nx,1,Ns,Nb]),[1 1 Nc 1 1]);
    im_currentseg_c=bsxfun(@times,im_currentseg,exp(1i*ph_c));
    ksp_currentseg_c=cfftn(im_currentseg_c,[1 2]);
    ksp_spen_c(m:Nseg:end,:,:,:,:,:)=reshape(A_currentseg*ksp_currentseg_c(:,:),dsize);
end

A_rec=crop(A,[Npe,Npe]);
AHA=A_rec'*A_rec;
[~,S,~]=svd(AHA);
lambda=0.5*S(1,1);
A_inv=pinv(AHA + lambda *eye(size(AHA)))*A_rec';
A_inv = pinv(A_rec);
A_1 = A(:,1:Npe); A_3 = A(:,2*Npe+1:end);
dsize=size(ksp_spen); 
switch add1
    case 1
        P = A_inv*A_1;
    case 2
        P = eye(Npe);
    case 3
        P = A_inv*A_3;
end
ksp_rec_c=reshape(P*A_inv*ksp_spen_c(:,:),dsize);
ksp_rec_c=bsxfun(@times,ksp_rec_c,exp(1i*2*pi*SPENpara.recon.PEshiftnorm*(-Npe/2:(Npe/2-1))).');
% im_rec_c=sos(flip(cifftn(ksp_rec_c,[1 2]),1),3);
disp('Doing SPEN reconstruction, done');
%%
ksize=[5 5];
im_res=zeros(Npe,Nx,Ns,Nb);
for m=1:Ns
    eigThresh_k = 0.05;
    eigThresh_im = 0.95;
    % crop a calibration area
    calib = crop(ksp_rec_c(:,:,:,m,1),[24,24,Nc]);
    [k,S] = dat2Kernel(calib,ksize);
    idx = max(find(S >= S(1)*eigThresh_k));
    [M,W] = kernelEig(k(:,:,:,1:idx),[Npe,Nx]);
    maps = M(:,:,:,end-1:end);
    % Weight the eigenvectors with soft-senses eigen-values
    weights = W(:,:,end-1:end) ;
    weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end-1:end) > eigThresh_im);
    weights = -cos(pi*weights)/2 + 1/2;
    % create and ESPIRiT operator
    ESP = ESPIRiT(maps,weights);
    for n=1:Nb
        nIterCG = 12; 
        disp('Performing ESPIRiT reconstruction from 2 maps')
        tic; [reskESPIRiT, resESPIRiT] = cgESPIRiT(ksp_rec_c(:,:,:,m,n),ESP, nIterCG, 0.01,ksp_rec_c(:,:,:,m,n)*0); toc
        im_res(:,:,m,n)=resESPIRiT(:,:,2);
    end
end
%%
% %% show the images
% im_final=round(2000*squeeze(im_rec_c)/max(im_rec_c(:)));
% disp('Showing the images after reconstruction');
% figure();imshowMRI(im_final,0,[Nb Ns]);
%%
im_final=round(2000*squeeze(flip(im_res,1))/max(im_res(:)));
winf=(hann(Npe)*hann(Nx).').^(1/3);
im_final_w=imresize(cifftn(bsxfun(@times,cfftn(im_final,[1 2]),winf),[1 2]),[Npe*4,Nx*4],'lanczos2');
disp('done.')
% disp('Showing the images after reconstruction');
% figure();imshowMRI(abs(im_final_w),0,[Nb Ns]);
