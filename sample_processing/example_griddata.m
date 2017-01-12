function []=example_griddata(choiceV3orV4);
%object:    a griddata example within gcmfaces
%inputs:    choiceV3orV4 ('v3' or 'v4') selects the sample GRID

gcmfaces_global;
if myenv.verbose>0;
    gcmfaces_msg('===============================================');
    gcmfaces_msg(['*** entering example_griddata that will define ' ...
        'a 2x2 gridded variable (an SSH field here), and then ' ...
        'use griddata to map it to the chosen grid. '],'');
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
original_Depth=mygrid.Depth;
original_Depth(original_Depth==0)=NaN;
lon=mygrid.XC; lat=mygrid.YC;
for ff=1:mygrid.nFaces;
    lon{ff}(1:2:end,:)=NaN;
    lat{ff}(1:2:end,:)=NaN;
    original_Depth{ff}(1:2:end,:)=NaN;
    lon{ff}(:,1:2:end)=NaN;
    lat{ff}(:,1:2:end)=NaN;
    original_Depth{ff}(:,1:2:end)=NaN;
end;
lon=convert2vector(lon);
lat=convert2vector(lat);
vec=convert2vector(original_Depth);

%%%%%%%%%%%%%%%%%%%%%%%
%make sure to accomodate the two longitude conventions:
if myenv.verbose>0;
    gcmfaces_msg('* accomodate both lon. conventions, then convert to data vector');
end;
x=[lon-360;lon;lon+360]; y=[lat;lat;lat]; z=[vec;vec;vec];
jj=find(~isnan(z)); x=x(jj); y=y(jj); z=z(jj);

%%%%%%%%%%%%%%%%%%%%%%%
%do the interpolation:
if myenv.verbose>0;
    gcmfaces_msg('* map/interpolate data using griddata');
end;
griddata_Depth=gcmfaces(5);
for ii=1:mygrid.nFaces;
    xi=mygrid.XC{ii}; yi=mygrid.YC{ii};
    zi = griddata(x',y',z',xi,yi);
    griddata_Depth{ii}=zi;
end;

%msk=mygrid.hFacC(:,:,1); griddata_Depth(find(msk==0))=NaN;

%%%%%%%%%%%%%%%%%%%%%%%
%illustrate the result:

figureL;
subplot(2,1,1); set(gca,'FontSize',14);
[X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,original_Depth);
pcolor(X,Y,FLD); axis([-180 180 -90 90]); shading flat;
title('input to griddata');
subplot(2,1,2); set(gca,'FontSize',14);
[X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,griddata_Depth);
pcolor(X,Y,FLD); axis([-180 180 -90 90]); shading flat;
title('output of griddata');

if myenv.verbose>0;
    gcmfaces_msg('* note of caution : griddata puts values on land. ');
    gcmfaces_msg('*** leaving example_griddata');
    gcmfaces_msg('===============================================');
end;
