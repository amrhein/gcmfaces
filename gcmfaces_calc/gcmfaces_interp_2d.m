function [val]=gcmfaces_interp_2d(fld,lon,lat,varargin);
%object:    linearly interpolate field to given positions
%inputs:    fld is a 2D gcmfaces field
%           lon,lat are position vectors
%optional:  doNearFill (1 by default) to use the nearest 
%               neighbor to complement linear interp
%outputs:   val is the vector of interpolated values
%
%pre-requisite: generate the delaunay triangulation using gcmfaces_bindata
%assumption: fld should show NaN for missing values

warning('off','MATLAB:dsearch:DeprecatedFunction');

global mygrid mytri myenv;

%inputs and pre-requisites
if ~isfield(myenv,'useDelaunayTri');
    myenv.useDelaunayTri=~isempty(which('DelaunayTri'));
end;

if isempty(mytri); 
    error('missing triangulation (mytri; from gcmfaces_bindata)');
end;

if nargin>3;
    doNearFill=varargin{1};
else;
    doNearFill=1;
end;

%switch longitude range to -180+180 or 0-360 according to grid
if max(mygrid.XC)<0;
    lon(find(lon>180))=lon(find(lon>180))-360;
end;

if max(mygrid.XC)>180;
    lon(find(lon<0))=lon(find(lon<0))+360;
end;

%do the actual interpolation    
if ~myenv.useDelaunayTri;%(code from old griddata.m)

    % Find the nearest triangle (t)
    x=convert2array(mygrid.XC); x=x(mytri.kk);
    y=convert2array(mygrid.YC); y=y(mytri.kk);
    VEC=convert2array(fld); VEC=VEC(mytri.kk);
    t = tsearch(x,y,mytri.TRI,lon',lat')';%the order of dims matters!!
    
    % Only keep the relevant triangles.
    out = find(isnan(t));
    if ~isempty(out), t(out) = ones(size(out)); end
    tri = mytri.TRI(t(:),:);
    
    % Compute Barycentric coordinates (w).  P. 78 in Watson.
    del = (x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1))) - ...
        (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1)));
    w(:,3) = ((x(tri(:,1))-lon(:)).*(y(tri(:,2))-lat(:)) - ...
        (x(tri(:,2))-lon(:)).*(y(tri(:,1))-lat(:))) ./ del;
    w(:,2) = ((x(tri(:,3))-lon(:)).*(y(tri(:,1))-lat(:)) - ...
        (x(tri(:,1))-lon(:)).*(y(tri(:,3))-lat(:))) ./ del;
    w(:,1) = ((x(tri(:,2))-lon(:)).*(y(tri(:,3))-lat(:)) - ...
        (x(tri(:,3))-lon(:)).*(y(tri(:,2))-lat(:))) ./ del;
    
    val = sum(VEC(tri) .* w,2);    
    val(out)=NaN;
    val=reshape(val,size(lon));
    
else;%(use TriScatteredInterp)
    
    %get interpolant
    VEC=convert2array(fld); VEC=VEC(mytri.kk);
    F = TriScatteredInterp(mytri.TRI, VEC);
    %do interpolation
    val=F(lon,lat);
    
end;

if doNearFill;%use the nearest neighbor to complement linear interp
    kk=gcmfaces_bindata(lon(:),lat(:));
    ARR=convert2array(fld);
    val2=ARR(kk);
    val(isnan(val))=val2(isnan(val));
end;


