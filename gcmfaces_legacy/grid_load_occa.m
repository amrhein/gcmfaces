

global mygrid;

if isempty(mygrid); 

dirGrid='/net/ross/raid2/gforget/1x1_50levels/GRID/';
list0={'XC','XG','YC','YG','RC','RF','RAC','DRC','DRF',...
        'DXC','DXG','DYC','DYG','hFacC','hFacS','hFacW','Depth'};
%list0={'AngleCS','AngleSN','Depth','DRC','DRF','DXC','DXG','DYC','DYG',...
%       'hFacC','hFacS','hFacW','maskCtrlC','maskCtrlS','maskCtrlW',...
%       'PHrefC','PHrefF','RAC','RAS','RAW','RAZ','RC','RF','XC','XG','YC','YG'};

list0={'XC','XG','YC','YG','RAC','DXC','DXG','DYC','DYG',...
       'hFacC','hFacS','hFacW','Depth'};
for iFld=1:length(list0);
   tmp1=gcmfaces(1,'ll'); eval(['tmp2=rdmds([dirGrid ''' list0{iFld} '*'']);']);
   tmp1.f1=tmp2; eval(['mygrid.' list0{iFld} '=tmp1;']);
end;

mygrid.AngleCS=mygrid.XC; mygrid.AngleCS(:,:)=1;
mygrid.AngleSN=mygrid.XC; mygrid.AngleSN(:,:)=0;

list0={'RC','RF','DRC','DRF'};
for iFld=1:length(list0);
   eval(['mygrid.' list0{iFld} '=rdmds([dirGrid ''' list0{iFld} '*'']);']);
end;

mygrid.hFacCsurf=mygrid.hFacC;
for ff=1:mygrid.hFacC.nFaces; mygrid.hFacCsurf{ff}=mygrid.hFacC{ff}(:,:,1); end;

end;

