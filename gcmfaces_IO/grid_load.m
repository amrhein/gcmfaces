function []=grid_load(dirGrid,nFaces,fileFormat,memoryLimit,omitNativeGrid);
%object:    load grid information, convert it to gcmfaces format
%           and encapsulate it in the global mygrid structure.
%inputs:    dirGrid is the directory where the grid files (gcm output) can be found.
%           nFaces is the number of faces in this gcm set-up of current interest.
%           fileFormat is the file format ('straight','cube','compact')
%optional:  memoryLimit is a flag that allows the user to omit secondary
%               grid fields in case memory/storage become an issue. It
%               takes 3 values : (0; the default) includes everything
%               (1) omits all 3D fields but hFacC (2) only loads XC & YC.
%           omitNativeGrid is a flag (0 by default) to bypass (when flag 1)
%               grid_load_native and grid_load_native_RAZ calls (that can
%               complement mygrid based upon e.g. tile*.mitgrid files)

if isempty(whos('memoryLimit')); memoryLimit=0; end;
if isempty(whos('omitNativeGrid')); omitNativeGrid=0; end;

gcmfaces_global; mygrid=[];

mygrid.dirGrid=dirGrid;
mygrid.nFaces=nFaces;
mygrid.fileFormat=fileFormat;
mygrid.gcm2facesFast=false;
mygrid.memoryLimit=memoryLimit;

if mygrid.memoryLimit>0;
    gcmfaces_msg(['* Warning from grid_load : memoryLimit>0 ' ...
        'may precludes advanced gmcfaces functions.'],'');
end;
if mygrid.memoryLimit>1;
    gcmfaces_msg(['* Warning from grid_load : memoryLimit>1 ' ...
        'may only allow for basic fields displays.'],'');
end;

if ~isempty(dir([dirGrid 'grid.specs.mat']));
    specs=open([dirGrid 'grid.specs.mat']);
    mygrid.ioSize=specs.ioSize;
    mygrid.facesSize=specs.facesSize;
    mygrid.facesExpand=specs.facesExpand;
    %example for creating grid.specs.mat, to put in dirGrid :
    %ioSize=[364500 1];
    %facesSize=[[270 450];[0 0];[270 270];[180 270];[450 270]];
    %facesExpand=[270 450];
    %save grid.specs.mat ioSize facesSize facesExpand;
elseif strcmp(fileFormat,'compact');
    v0=rdmds([dirGrid 'XC']);
    mygrid.ioSize=size(v0);
    nn=size(v0,1); pp=size(v0,2)/nn;
    mm=(pp+4-mygrid.nFaces)/4*nn;
    mygrid.facesSize=[[nn mm];[nn mm];[nn nn];[mm nn];[mm nn];[nn nn]];
    mygrid.facesExpand=[];
elseif strcmp(fileFormat,'cube');
    v0=rdmds([dirGrid 'XC']);
    mygrid.ioSize=size(v0);
    nn=size(v0,2);
    mygrid.facesSize=[[nn nn];[nn nn];[nn nn];[nn nn];[nn nn];[nn nn]];
    mygrid.facesExpand=[];
elseif strcmp(fileFormat,'straight');
    v0=rdmds([dirGrid 'XC']);
    mygrid.ioSize=size(v0);
    mygrid.facesSize=mygrid.ioSize;
    mygrid.facesExpand=[];
end;
mygrid.missVal=NaN;%will be set to 0 once the grid has been loaded.

if  ~(nFaces==1&strcmp(fileFormat,'straight'))&...
        ~(nFaces==6&strcmp(fileFormat,'cube'))&...
        ~(nFaces==6&strcmp(fileFormat,'compact'))&...
        ~(nFaces==5&strcmp(fileFormat,'compact'));
    %     if myenv.verbose;
    %         fprintf('\nconvert2gcmfaces.m init: there are several supported file conventions. \n');
    %         fprintf('  By default gcmfaces assumes MITgcm type binary formats as follows: \n')
    %         fprintf('  (1 face) straight global format; (4 or 5 faces) compact global format\n');
    %         fprintf('  (6 faces) cube format with one face after the other. \n');
    %         fprintf('  If this is inadequate, you can change the format below.\n\n');
    %     end;
    error('non-tested topology/fileFormat');
end;

%the various grid fields:
if mygrid.memoryLimit==0;
    list0={'XC','XG','YC','YG','RAC','RAZ','DXC','DXG','DYC','DYG','hFacC','hFacS','hFacW','Depth'};
elseif mygrid.memoryLimit==1;
    list0={'XC','XG','YC','YG','RAC','RAZ','DXC','DXG','DYC','DYG','hFacC','Depth'};
elseif mygrid.memoryLimit==2;
    list0={'XC','YC'};
end;

for iFld=1:length(list0);
    eval(['mygrid.' list0{iFld} '=rdmds2gcmfaces([dirGrid ''' list0{iFld} '*'']);']);
end;

%the vertical grid
list0={'RC','RF','DRC','DRF'};
for iFld=1:length(list0);
    eval(['mygrid.' list0{iFld} '=squeeze(rdmds([dirGrid ''' list0{iFld} '*'']));']);
end;

%grid orientation
if mygrid.memoryLimit<2;
    list0={'AngleCS','AngleSN'};
    test0=~isempty(dir([dirGrid 'AngleCS*']));
    if test0;
        for iFld=1:length(list0);
            eval(['mygrid.' list0{iFld} '=rdmds2gcmfaces([dirGrid ''' list0{iFld} '*'']);']);
        end;
    else;
        warning('\n AngleCS/AngleSN not found; set to 1/0 assuming lat/lon grid.\n');
        mygrid.AngleCS=mygrid.XC; mygrid.AngleCS(:)=1;
        mygrid.AngleSN=mygrid.XC; mygrid.AngleSN(:)=0;
    end;
end;

%if grid is incomplete (blank tiles) then try to get 
%additional info from native grid, or apply missing value.
test0=sum(isnan(mygrid.XC))>0;
test1=prod(mygrid.ioSize)~=sum(isnan(NaN*mygrid.XC(:)));
if (test0|test1)&~omitNativeGrid;
  %treat fields that are part of the native grid
  mygrid1=mygrid; mygrid=[];
  grid_load_native(dirGrid,nFaces,0);
  mygrid1.XC=mygrid.XC; mygrid1.YC=mygrid.YC;
  mygrid1.XG=mygrid.XG; mygrid1.YG=mygrid.YG;
  mygrid1.RAC=mygrid.RAC;
  mygrid=mygrid1;
  %apply missing value for fields that aren't
  list0={'hFacC','hFacS','hFacW','Depth','AngleCS','AngleSN'};
  for ii=1:length(list0); 
    eval(['tmp1=mygrid.' list0{ii} ';']);
    tmp1(isnan(tmp1))=0;
    eval(['mygrid.' list0{ii} '=tmp1;']);
  end;
  %and fix angles if needed
  tmp1=mygrid.AngleCS.^2+mygrid.AngleSN.^2;
  tmp1=1*(tmp1>0.999&tmp1<1.001);
  mygrid.AngleCS(tmp1==0)=1;
  mygrid.AngleSN(tmp1==0)=0;
end;

%get full RAZ (incl. 'extra line and column') needed for e.g. rotational computations
if mygrid.memoryLimit<2&~omitNativeGrid;
    grid_load_native_RAZ;
end;

%grid masks
if mygrid.memoryLimit<1;
    mygrid.hFacCsurf=mygrid.hFacC;
    for ff=1:mygrid.hFacC.nFaces; mygrid.hFacCsurf{ff}=mygrid.hFacC{ff}(:,:,1); end;
    
    mskC=mygrid.hFacC; mskC(mskC==0)=NaN; mskC(mskC>0)=1; mygrid.mskC=mskC;
    mskW=mygrid.hFacW; mskW(mskW==0)=NaN; mskW(mskW>0)=1; mygrid.mskW=mskW;
    mskS=mygrid.hFacS; mskS(mskS==0)=NaN; mskS(mskS>0)=1; mygrid.mskS=mskS;
end;

%zonal mean and sections needed for transport computations
if mygrid.memoryLimit<1;
    if ~isfield(mygrid,'mygrid.LATS_MASKS');
        gcmfaces_lines_zonal;
        mygrid.LATS=[mygrid.LATS_MASKS.lat]';
    end;
    if ~isfield(mygrid,'LINES_MASKS');
        [lonPairs,latPairs,names]=line_greatC_TUV_MASKS_v4;
        gcmfaces_lines_transp(lonPairs,latPairs,names);
    end;
end;

%to allow convert2gcmfaces/doFast:
if isempty(mygrid.facesExpand)&mygrid.memoryLimit<2;
    tmp1=convert2gcmfaces(mygrid.XC);
    tmp1(:)=[1:length(tmp1(:))];
    nn=length(tmp1(:));
    mygrid.gcm2faces=convert2gcmfaces(tmp1);
    mygrid.faces2gcmSize=size(tmp1);
    mygrid.faces2gcm=convert2gcmfaces(tmp1);
    for iFace=1:mygrid.nFaces;
        n=length(mygrid.gcm2faces{iFace}(:));
        mygrid.faces2gcm{iFace}=mygrid.gcm2faces{iFace}(:);
        mygrid.gcm2faces{iFace}=sparse([1:n],mygrid.gcm2faces{iFace}(:),ones(1,n),n,nn);
    end;
    mygrid.gcm2facesFast=true;
end;

%reset missVal parameter to 0. 
%Note : this is only used by convert2widefaces, for runs with cropped grids.
%Note : 0 should not be used as a fill for the grid itself (NaN was used).
mygrid.missVal=0;

