function []=example_bin_average(choiceV3orV4,doSmooth);
%object:    a bin average example within gcmfaces
%inputs:    choiceV3orV4 ('v3' or 'v4') selects the sample GRID
%           doSmooth; if set to 1, follow on with smoothing example

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

mygrid=[]; grid_load(dirGrid,nF,fileFormat,0);

if myenv.verbose>0;
    gcmfaces_msg('===============================================');
%     gcmfaces_msg('*** entering example_bin_average that will define','');
%     gcmfaces_msg('randomly distributed variable (position and value),','  ');
%     gcmfaces_msg('then bin average it to the chosen grid,','  ');
%     gcmfaces_msg('and apply a smoothing filter to it','  ');
    gcmfaces_msg(['*** entering example_bin_average that will first define ' ... 
        'randomly distributed variable (position and value), ' ... 
        'then bin average it to the chosen grid, ' ... 
        'and finally apply a smoothing filter to it '],'');
end;

warning('off','MATLAB:dsearch:DeprecatedFunction');
warning('off','MATLAB:delaunay:DuplicateDataPoints');

%%%%%%%%%%%%%%%%%%%%%%%
%generate random data

if myenv.verbose>0;
    gcmfaces_msg('* generate random data');
end;
nobs=1e6;
lat=(rand(nobs,1)-0.5)*2*90;
lon=(rand(nobs,1)-0.5)*2*180; 
%needed for 0-360 longitude convention
if mygrid.nFaces==1;
    xx=find(lon<0);lon(xx)=lon(xx)+360;
end;
obsPts=(rand(nobs,1)-0.5)*2;

%%%%%%%%%%%%%%%%%%%%%%%
%generate delaunay triangulation

if myenv.verbose>0;
    gcmfaces_msg('* call gcmfaces_bindata : generate delaunay triangulation');
end;
gcmfaces_bindata;

%%%%%%%%%%%%%%%%%%%%%%%%
%bin average random data

if myenv.verbose>0;
    gcmfaces_msg('* call gcmfaces_bindata : bin average data');
end;
obsMap=gcmfaces_bindata(lon,lat,obsPts);

%%%%%%%%%%%%%%%%%%%%%%%%
%format conversions
if myenv.verbose>0;
    gcmfaces_msg('* call convert2array : convert from gcmfaces to array format');
end;
obsArray=convert2array(obsMap);%put in gcmfaces format
if myenv.verbose>0;
    gcmfaces_msg('* call convert2array : convert back to gcmfaces');
end;
obsMap2=convert2array(obsArray);%put in gcmfaces format
if myenv.verbose>0;
    gcmfaces_msg('* call convert2gcmfaces : convert from gcmfaces to file format');
end;
obsOut=convert2gcmfaces(obsMap);	%put in gcm input format
if myenv.verbose>0;
    gcmfaces_msg('* summarizing data formats:');
    aa=whos('obs*');
    aaa=aa(4); bb=['[' num2str(aaa.size(1)) 'x' num2str(aaa.size(2)) ']'];
    bb=fprintf('   %8s %8s %12s    : data points (vector)\n',aaa.name,aaa.class,bb);
    aaa=aa(2); bb=['[' num2str(aaa.size(1)) 'x' num2str(aaa.size(2)) ']'];
    bb=fprintf('   %8s %8s %12s    : gridded data (gcmfaces)\n',aaa.name,aaa.class,bb);
    aaa=aa(1); bb=['[' num2str(aaa.size(1)) 'x' num2str(aaa.size(2)) ']'];
    bb=fprintf('   %8s %8s %12s    : array format\n',aaa.name,aaa.class,bb);
    aaa=aa(3); bb=['[' num2str(aaa.size(1)) 'x' num2str(aaa.size(2)) ']'];
    bb=fprintf('   %8s %8s %12s    : output format\n',aaa.name,aaa.class,bb);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%quick display
if myenv.verbose>0;
    gcmfaces_msg('* crude display of results in array format');
end;
figureL; imagescnan(obsArray','nancolor',[1 1 1]*0.8); axis xy; caxis([-1 1]*0.4); colorbar;
title('bin averaged data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now illustrate the smoothing filter
if myenv.verbose>0;
    gcmfaces_msg('* entering example_smooth : apply smoother to gridded data');
end;
if doSmooth; example_smooth; end;
if myenv.verbose>0;
    gcmfaces_msg('* leaving example_smooth');
end;

if myenv.verbose>0;
    gcmfaces_msg('*** leaving example_bin_average');
    gcmfaces_msg('===============================================');
end;




