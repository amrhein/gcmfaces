function []=example_faces2latlon2faces(choiceV3orV4);
%object:    interpolate back and forth between gcmfaces and lat-lon grid

if 0;
    clear all; grid_load('./',6,'cube'); gcmfaces_bindata; gcmfaces_global;
    lon=[-179:2:179]; lat=[-30:2:30];
    [lat,lon] = meshgrid(lat,lon);
    fld=sin(2*pi/180*mygrid.XC).*sin(2*pi/180*mygrid.YC);
    fld2=gcmfaces_interp_2d(fld,lon,lat);%go from gcmfaces grid to lat-lon grid
end;

gcmfaces_global; global mytri;

if myenv.verbose>0;
    gcmfaces_msg('===============================================');
    gcmfaces_msg(['*** entering example_faces2latlon2faces that will ' ...
        'interpolate a field on the chosen to a lat-lon grid and back' ...
        'to the original grid'],'');
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

%store reference grid
mygrid_refgrid=mygrid;

%%%%%%%%%%%%%%%%%%%
%field for testing:
%%%%%%%%%%%%%%%%%%%

if myenv.verbose>0;
    gcmfaces_msg('* define a sinusoidal field for testing');
end;
fld=sin(2*pi/180*mygrid.XC).*sin(2*pi/180*mygrid.YC);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%lat-lon grid for testing:
%%%%%%%%%%%%%%%%%%%%%%%%%%

if myenv.verbose>0;
    gcmfaces_msg('* define 2x2 lat-lon grid for testing');
end;
%define lat-lon grid
lon=[-179:2:179]; lat=[-89:2:89]; aa=[-180 180 -90 90];
if max(mygrid_refgrid.XC)>180;
    lon=[1:2:359]; aa=[0 360 -90 90];
end;
% lon=[-179:2:179]; lat=[-30:2:30];
[lat,lon] = meshgrid(lat,lon);
%prepare mygrid for lat-lon with no mask
mygrid_latlon.nFaces=1;
mygrid_latlon.XC=gcmfaces({lon}); mygrid_latlon.YC=gcmfaces({lat});
mygrid_latlon.dirGrid='none';
mygrid_latlon.fileFormat='straight';
mygrid_latlon.ioSize=size(lon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interpolate to lat-lon grid:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if myenv.verbose>0;
    gcmfaces_msg('* interpolate test field to lat-lon grid');
end;
mygrid=mygrid_latlon; gcmfaces_bindata; 
veclon=convert2array(mygrid.XC); veclon=veclon(mytri.kk);
veclat=convert2array(mygrid.YC); veclat=veclat(mytri.kk);

mygrid=mygrid_refgrid; gcmfaces_bindata; 
vecfld=gcmfaces_interp_2d(fld,veclon,veclat);%go from gcmfaces grid to lat-lon grid

mygrid=mygrid_latlon; gcmfaces_bindata; 
fld_latlon=NaN*convert2array(mygrid.XC); 
fld_latlon(mytri.kk)=vecfld;
fld_latlon=convert2array(fld_latlon); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interpolate back to ref grid:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if myenv.verbose>0;
    gcmfaces_msg('* interpolate back to original grid');
end;
mygrid=mygrid_refgrid; gcmfaces_bindata; 
veclon=convert2array(mygrid.XC); veclon=veclon(mytri.kk);
veclat=convert2array(mygrid.YC); veclat=veclat(mytri.kk);

mygrid=mygrid_latlon; gcmfaces_bindata; 
vecfld=gcmfaces_interp_2d(fld_latlon,veclon,veclat);%go from gcmfaces grid to lat-lon grid

mygrid=mygrid_refgrid; gcmfaces_bindata; 
fld_refgrid=NaN*convert2array(mygrid.XC); 
fld_refgrid(mytri.kk)=vecfld;
fld_refgrid=convert2array(fld_refgrid); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get and plot the lat-lon array:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if myenv.verbose>0;
    gcmfaces_msg('* plot lat-lon grid interpolated field');
end;
FLD=fld_latlon{1};
figureL; 
subplot(2,1,1); set(gca,'FontSize',14);
pcolor(lon,lat,FLD); axis([-270 270 -100 100]); 
shading flat; colorbar; title('interpolated lat-lon map');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot difference due to interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if myenv.verbose>0;
    gcmfaces_msg('* plot error due to back and forth interpolation');
end;
mygrid=mygrid_refgrid;
subplot(2,1,2); set(gca,'FontSize',14);
[X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,fld_refgrid-fld);
pcolor(X,Y,log10(abs(FLD))); axis([-270 270 -100 100]); 
shading flat; caxis([-5 -2]); colorbar; 
title('log10(error) due to back and forth interpolation');

if myenv.verbose>0;
    gcmfaces_msg('*** leaving example_faces2latlon2faces');
    gcmfaces_msg('===============================================');
end;
