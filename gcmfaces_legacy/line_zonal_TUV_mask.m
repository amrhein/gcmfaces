function [line_out]=line_zonal_TUV_mask(lat_in);

global mygrid;

y=mygrid.YC;
yP1=exch_T_N(y);

%T part:
yy=lat_in;
mm=yP1<=yy; mm(isnan(mm))=0;
for iF=1:y.nFaces;
   eval(['tmp1=mm.f' num2str(iF) ';']); 
   tmp2=circshift(tmp1,[-1 0])+circshift(tmp1,[1 0])+...
        circshift(tmp1,[0 -1])+circshift(tmp1,[0 1]);
   tmp2(tmp1>0)=0; tmp2(tmp2~=0)=1; tmp2(tmp2==0)=NaN;
   eval(['mm.f' num2str(iF) '=tmp2(2:end-1,2:end-1);']);
end;
mmt=mm;
%UV part:
nn=exch_T_N(mmt);
mm=yP1<=yy; mm(isnan(mm))=0;
mmu=y; mmv=y; 
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
%store:
line_out=struct('mmt',mmt,'mmu',mmu,'mmv',mmv,'lat',lat_in);


