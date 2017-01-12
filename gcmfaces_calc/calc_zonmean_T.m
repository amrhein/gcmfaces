function [FLD,X,Y]=calc_zonmean_T(fld);
%object:    compute zonal mean
%inputs:    fld is the field of interest
%output:    FLD is the zonal mean field
%
%notes:     mygrid.LATS_MASKS is the set of quasi longitudinal lines along which
%               means will computed, as computed in gcmfaces_lines_zonal

global mygrid;

%initialize output:
n3=max(size(fld.f1,3),1); n4=max(size(fld.f1,4),1);
FLD=NaN*squeeze(zeros(length(mygrid.LATS_MASKS),n3,n4));

%use array format to speed up computation below:
fld=convert2array(fld);
n1=size(fld,1); n2=size(fld,2);
fld=reshape(fld,n1*n2,n3*n4);
rac=reshape(convert2array(mygrid.RAC),n1*n2,1)*ones(1,n3*n4);
if n3==length(mygrid.RC);
    hFacC=reshape(convert2array(mygrid.hFacC),n1*n2,n3*n4);
else;
    hFacC=reshape(convert2array(mygrid.mskC(:,:,1)),n1*n2,1)*ones(1,n3*n4);
    hFacC(isnan(hFacC))=0;
end;
%masked area only:
rac(isnan(fld))=0;
hFacC(isnan(fld))=0;
fld(isnan(fld))=0;

rac=rac.*hFacC;
fld=fld.*rac;

for iy=1:length(mygrid.LATS_MASKS);
    
    %get list ofpoints that form a zonal band:
    mm=convert2array(mygrid.LATS_MASKS(iy).mskCedge);
    mm=find(~isnan(mm)&mm~=0);
    
    %do the area weighed average along this band:
    tmp1=sum(fld(mm,:),1);
    tmp2=sum(rac(mm,:),1);
    tmp2(tmp2==0)=NaN;
    tmp1=tmp1./tmp2;
    
    %store:
    FLD(iy,:,:)=reshape(tmp1,n3,n4);
    
end;

X=[]; Y=[];
if size(FLD,2)==length(mygrid.RC);
    X=mygrid.LATS*ones(1,length(mygrid.RC));
    Y=ones(length(mygrid.LATS),1)*(mygrid.RC');
elseif size(FLD,2)==1;
    X=mygrid.LATS;
    Y=ones(length(mygrid.LATS),1);
end;

