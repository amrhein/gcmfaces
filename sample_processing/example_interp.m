function []=example_interp(choiceV3orV4);
%object:    a interp2 example within gcmfaces
%inputs:    choiceV3orV4 ('v3' or 'v4') selects the sample GRID

gcmfaces_global;
if myenv.verbose>0;
    gcmfaces_msg(['* call example_interp : will define ' ...
        'a 2x2 gridded variable (an SSH field here), and then ' ...
        'use interp2 to map it to the chosen grid. '],'');
end;

%%%%%%%%%%%%%%%%%
%load parameters:
%%%%%%%%%%%%%%%%%

gcmfaces_global;
dir0=[myenv.gcmfaces_dir '/sample_input/'];
if strcmp(choiceV3orV4,'v4');
    nF=5;
    fileFormat='compact';
    dirGrid=[dir0 '../../GRID/'];
else;
    error('files are missing');
    nF=1;
    fileFormat='straight';
    dirGrid=[dir0 '/GRID' choiceV3orV4  '/'];
end;

mygrid=[]; grid_load(dirGrid,nF,fileFormat,1);

%%%%%%%%%%%%%%%%%%%%%%%
%get sample data: subsampled mygrid.Depth
if myenv.verbose>0;
    gcmfaces_msg('* define test case : subsampled bathymetry');
end;
lon=mygrid.XC; lat=mygrid.YC;
original_Depth=mygrid.Depth;
original_Depth(original_Depth==0)=NaN;
%
[lon,lat,original_Depth]=convert2pcol(lon,lat,original_Depth);
nn=[62:239]; lon=lon(:,nn); lat=lat(:,nn);
original_Depth=original_Depth(:,nn);
[lon,lat]=meshgrid(lon(:,1),lat(1,:));
original_Depth=original_Depth';

%%%%%%%%%%%%%%%%%%%%%%%
%do the interpolation:

interp2_depth=gcmfaces(5);
for ii=1:mygrid.nFaces;
    xi=mygrid.XC{ii}; yi=mygrid.YC{ii};
    zi = interp2(lon,lat,original_Depth,xi,yi);
    interp2_depth{ii}=zi;
end;

%%%%%%%%%%%%%%%%%%%%%%%
%illustrate the result:

figureL;
subplot(2,1,1); set(gca,'FontSize',14);
pcolor(lon,lat,original_Depth); axis([-180 180 -90 90]); shading flat;
title('input to interp2');
subplot(2,1,2); set(gca,'FontSize',14);
[X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,interp2_depth);
pcolor(X,Y,FLD); axis([-180 180 -90 90]); shading flat;
title('output of interp2');


