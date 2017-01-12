function [line_out]=line_greatC_TUV_mask(lonPair_in,latPair_in);

global mygrid;

%get carthesian coordinates:
%... of grid
lon=mygrid.XC; lat=mygrid.YC;
x=cos(lat*pi/180).*cos(lon*pi/180);
y=cos(lat*pi/180).*sin(lon*pi/180); 
z=sin(lat*pi/180); 
%... and of end points 
x0=cos(latPair_in*pi/180).*cos(lonPair_in*pi/180); 
y0=cos(latPair_in*pi/180).*sin(lonPair_in*pi/180);
z0=sin(latPair_in*pi/180);                 

%get the rotation matrix:
%1) rotate around x axis to put first point at z=0 
theta=atan2(-z0(1),y0(1));
R1=[[1;0;0] [0;cos(theta);sin(theta)] [0;-sin(theta);cos(theta)]];
tmp0=[x0;y0;z0]; tmp1=R1*tmp0; x1=tmp1(1,:); y1=tmp1(2,:); z1=tmp1(3,:);
x0=x1; y0=y1; z0=z1; 
%2) rotate around z axis to put first point at y=0 
theta=atan2(x0(1),y0(1));
R2=[[cos(theta);sin(theta);0] [-sin(theta);cos(theta);0] [0;0;1]];
tmp0=[x0;y0;z0]; tmp1=R2*tmp0; x1=tmp1(1,:); y1=tmp1(2,:); z1=tmp1(3,:);
x0=x1; y0=y1; z0=z1;
%3) rotate around y axis to put second point at z=0 
theta=atan2(-z0(2),-x0(2));
R3=[[cos(theta);0;-sin(theta)] [0;1;0] [sin(theta);0;cos(theta)]];
tmp0=[x0;y0;z0]; tmp1=R3*tmp0; x1=tmp1(1,:); y1=tmp1(2,:); z1=tmp1(3,:);
x0=x1; y0=y1; z0=z1;

%apply rotation to grid:
tmpx=convert2array(x); tmpy=convert2array(y); tmpz=convert2array(z);
tmp1=find(~isnan(tmpx));
tmpx2=tmpx(tmp1); tmpy2=tmpy(tmp1); tmpz2=tmpz(tmp1);
tmp2=[tmpx2';tmpy2';tmpz2'];
tmp3=R3*R2*R1*tmp2;
tmpx2=tmp3(1,:); tmpy2=tmp3(2,:); tmpz2=tmp3(3,:);
tmpx(tmp1)=tmpx2; tmpy(tmp1)=tmpy2; tmpz(tmp1)=tmpz2;
x=convert2array(tmpx); y=convert2array(tmpy); z=convert2array(tmpz);

%compute the great circle mask:
zz=exch_T_N(z);
mm=zz<=0; mm(isnan(mm))=0;
for iF=1:y.nFaces;
   eval(['tmp1=mm.f' num2str(iF) ';']); 
   tmp2=circshift(tmp1,[-1 0])+circshift(tmp1,[1 0])+...
        circshift(tmp1,[0 -1])+circshift(tmp1,[0 1]);
   tmp2(tmp1>0)=0; tmp2(tmp2~=0)=1; tmp2(tmp2==0)=NaN;
   eval(['mm.f' num2str(iF) '=tmp2(2:end-1,2:end-1);']);
end;

%UV part:
mmt=mm;
nn=exch_T_N(mmt);
mm=zz<=0; mm(isnan(mm))=0;
mmu=mmt; mmv=mmt; 
for iFace=1:y.nFaces; 
   iF=num2str(iFace);
   eval(['tmp_u=mm.f' iF '(1:end-1,:)==1&nn.f' iF '(2:end,:)==1;']);
   eval(['tmp_u2=mm.f' iF '(2:end,:)==1&nn.f' iF '(1:end-1,:)==1;']);
   tmp_u=1*(tmp_u(1:end-1,2:end-1)-tmp_u2(1:end-1,2:end-1)); eval(['mmu.f' iF '=tmp_u;']);
   eval(['tmp_v=mm.f' iF '(:,1:end-1)==1&nn.f' iF '(:,2:end)==1;']);
   eval(['tmp_v2=mm.f' iF '(:,2:end)==1&nn.f' iF '(:,1:end-1)==1;']);
   tmp_v=1*(tmp_v(2:end-1,1:end-1)-tmp_v2(2:end-1,1:end-1)); eval(['mmv.f' iF '=tmp_v;']);
end;
mmu(mmu==0)=NaN; mmv(mmv==0)=NaN;

for kk=1:3;
%select field to treat:
switch kk;
case 1; mm=mmt;
case 2; mm=mmu;
case 3; mm=mmv;
end;
%split in two contours: 
theta=[];
theta(1)=atan2(y0(1),x0(1));
theta(2)=atan2(y0(2),x0(2));

tmpx=convert2array(x); tmpy=convert2array(y); tmpz=convert2array(z);
tmptheta=atan2(tmpy,tmpx);
tmpm=convert2array(mm); tmpmIn=tmpm; tmpmOut=tmpmIn;
if theta(2)<0; 
   tmp00=find(tmptheta<=theta(2)); tmptheta(tmp00)=tmptheta(tmp00)+2*pi; 
   theta(2)=theta(2)+2*pi; 
end;
tmpmIn(find(tmptheta>theta(2)|tmptheta<theta(1)))=NaN;
tmpmOut(find(tmptheta<=theta(2)&tmptheta>=theta(1)))=NaN;
mmIn=convert2array(tmpmIn); mmOut=convert2array(tmpmOut);
%ensure that we select the shorther segment as mmIn:
if theta(2)-theta(1)>pi; tmp1=mmIn; mmIn=mmOut; mmOut=tmp1; end;
%store result:
switch kk;
case 1; mmtIn=mmIn; mmtOut=mmOut;
case 2; mmuIn=mmIn; mmuOut=mmOut;
case 3; mmvIn=mmIn; mmvOut=mmOut;
end;
end;

%output:
line_out=struct('mmt',mmt,'mmtIn',mmtIn,'mmtOut',mmtOut,...
   'mmu',mmu,'mmuIn',mmuIn,'mmuOut',mmuOut,...
   'mmv',mmv,'mmvIn',mmvIn,'mmvOut',mmvOut,...
   'lonPair',lonPair_in,'latPair',latPair_in);



