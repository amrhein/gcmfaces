

global mygrid;

if isempty(mygrid); 

dirGrid='/net/altix3700/raid4/gforget/mysetups/ecco_v4/RUNS/GRIDmds_90x50/';
list0={'XC','XG','YC','YG','RC','RF','RAC','DRC','DRF',...
        'DXC','DXG','DYC','DYG','hFacC','hFacS','hFacW','Depth'};
%list0={'AngleCS','AngleSN','Depth','DRC','DRF','DXC','DXG','DYC','DYG',...
%       'hFacC','hFacS','hFacW','maskCtrlC','maskCtrlS','maskCtrlW',...
%       'PHrefC','PHrefF','RAC','RAS','RAW','RAZ','RC','RF','XC','XG','YC','YG'};

list0={'XC','XG','YC','YG','RAC','DXC','DXG','DYC','DYG',...
       'hFacC','hFacS','hFacW','Depth','AngleCS','AngleSN'};
for iFld=1:length(list0);
   eval(['mygrid.' list0{iFld} '=rdmds_compact2llc([dirGrid ''' list0{iFld} '*'']);']);
end;

list0={'RC','RF','DRC','DRF'};
for iFld=1:length(list0);
   eval(['mygrid.' list0{iFld} '=rdmds([dirGrid ''' list0{iFld} '*'']);']);
end;

mygrid.hFacCsurf=mygrid.hFacC;
for ff=1:mygrid.hFacC.nFaces; mygrid.hFacCsurf{ff}=mygrid.hFacC{ff}(:,:,1); end;

end;

