function []=gcmfaces_lines_zonal(varargin);
%object:    define the set of quasi longitudinal lines
%           along which transports will integrated;
%           LATS_MASKS that will be added to mygrid.
%optional input:
%           LATS_VAL is the latitudes vector ([-89:89]' by default)
%           MSKS_NAM is the list of fields to store in mygrid

global mygrid;

if nargin>0; LATS_VAL=varargin{1}; else; LATS_VAL=[-89:89]'; end;
if nargin>1; MSKS_NAM=varargin{2}; else; MSKS_NAM={'mskCint','mskCedge','mskWedge','mskSedge'}; end;

for iy=1:length(LATS_VAL);

    mskCint=1*(mygrid.YC>=LATS_VAL(iy));
    [mskCedge,mskWedge,mskSedge]=gcmfaces_edge_mask(mskCint);

    for im=1:length(MSKS_NAM);
      eval(['tmp1.' MSKS_NAM{im} '=' MSKS_NAM{im} ';']);
    end;
    tmp1.lat=LATS_VAL(iy);

    %store:
    if iy==1;
        LATS_MASKS=tmp1;
    else;
        LATS_MASKS(iy)=tmp1;
    end;
    
end;

mygrid.LATS_MASKS=LATS_MASKS;


