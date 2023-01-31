function ksp=SPEN_parse_siemenstwix(twix_obj)

if iscell(twix_obj)
    twix_obj=twix_obj{end};
end

Nx=twix_obj.image.NCol/2;
Nc=twix_obj.image.NCha;
Netl=twix_obj.image.NLin;
% Nseg=twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{9};
Nseg = 1;   % zsj 0702
Ns=twix_obj.image.NSli;
if (isfield(twix_obj,'phasecor'))
    ksp_PC=reshape(twix_obj.phasecor.unsorted(),Nx,Nc,3,[]);
    ksp_PC=ksp_PC(:,:,[2 1 3],:);
    ksp_unsort=reshape(twix_obj.image.unsorted(),Nx,Nc,Netl,[]); 
    ksp_unsort_c=EPIEvenoddFix(ksp_unsort,ksp_PC,1,3);
end


ksp=reshape(ksp_unsort_c,Nx,Nc,Netl,Ns,Nseg,[]);
slicePos=twix_obj.image.slicePos;
[~,sliceorder]=sort(unique(slicePos(3,:),'stable'));
IndexCell=repmat({':'},[1 length(size(ksp))]);
IndexCell{4}=sliceorder;
ksp=ksp(IndexCell{:}); 

ncc=8;
porder=[1 3 2];
ksp_scc=reshape(ksp,Nx,Nc,[]);
ksp_scc=permute(ksp_scc,porder);
sccmtx = calcSCCMtx(ksp_scc);
ksp_scc = CC(ksp_scc,sccmtx(:,1:ncc));
ksp_scc=ipermute(ksp_scc,porder);
ksp=reshape(ksp_scc,Nx,ncc,Netl,Ns,Nseg,[]);
