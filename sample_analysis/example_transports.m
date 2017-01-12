function [diags]=example_transports(choiceV3orV4);
%object:    illustrates a series of standard computations
%           (streamfunctions, transports, zonal means, etc.)
%inputs:    choiceV3orV4 ('v3' or 'v4') selects the sample GRID

gcmfaces_global;
if myenv.verbose>0;
    gcmfaces_msg('===============================================');
    gcmfaces_msg(['*** entering example_transports ' ...
        'that will load a set of variables form file ' ...
        'and compute a series of derived diagnostics '],'');
    %     gcmfaces_msg('*** entering example_transports','');
    %     gcmfaces_msg('that will load a set of variables form file','  ');
    %     gcmfaces_msg('and compute a series of derived diagnostics','  ');
end;

%%%%%%%%%%%%%%%%%
%load parameters:
%%%%%%%%%%%%%%%%%

if myenv.verbose>0;
    gcmfaces_msg('* set grid files path and format, and number of faces for this grid');
end;
myenv.nctiles=1;
if strcmp(choiceV3orV4,'v4');
    nF=5;
    fileFormat='compact';
    dir0=fullfile([myenv.gcmfaces_dir '../']);
    dirGrid=fullfile(dir0,'GRID/');
    myenv.nctilesdir=fullfile(dir0,'release1/nctiles_climatology/');
else;
    error('files are missing');
    nF=1;
    dir0=fullfile(myenv.gcmfaces_dir,'sample_input/');
    fileFormat='straight';
    dirGrid=fullfile(dir0,'grid_occa/');
    myenv.nctilesdir=fullfile(dir0,'nctiles_occa/');
end;
if myenv.verbose>0;
    gcmfaces_msg('* call grid_load : load grid to memory (mygrid) according to');
    gcmfaces_msg(['dirGrid = ' dirGrid],'  ');
    gcmfaces_msg(['nFaces = ' num2str(nF)],'  ');
    gcmfaces_msg(['fileFormat = ' fileFormat],'  ');
    %     fprintf(['  > dirGrid = ' dirGrid '\n']);
end;

if ~isdir(myenv.nctilesdir);
    diags=[];
    warning(['skipping example_transports (missing ' myenv.nctilesdir ')']);
    return;
end;

% warning('skipping grid_load\n');
grid_load(dirGrid,nF,fileFormat);

if myenv.verbose>0;
    gcmfaces_msg('* call gcmfaces_lines_zonal : determine grid lines that closely follow');
    gcmfaces_msg('parallel lines and will be used in zonal mean and overturning computations','  ');
end;
% warning('skipping gcmfaces_lines_zonal\n');
if strcmp(choiceV3orV4,'v4');
    gcmfaces_lines_zonal;
else;
    gcmfaces_lines_zonal([-75:75]');
end;
if myenv.verbose>0;
    gcmfaces_msg('* call gcmfaces_lines_transp : determine grid lines that closely follow');
    gcmfaces_msg('great circles and will be used to compute transsects transports','  ');
end;
% warning('skipping gcmfaces_lines_transp\n');
eval(['[lonPairs,latPairs,names]=line_greatC_TUV_MASKS_' choiceV3orV4 ';']);
gcmfaces_lines_transp(lonPairs,latPairs,names);

%%%%%%%%%%%%%%%%%
%do computations:
%%%%%%%%%%%%%%%%%

% [listTimes]=diags_list_times;
diags.listTimes=1;

if myenv.verbose>0; gcmfaces_msg('* call rdmds2gcmfaces : load velocity fields');end;
listVars={'UVELMASS','VVELMASS'};
listDiags={'fldBAR','gloOV','fldTRANSPORTS','gloMT_FW'};
for vvv=1:length(listVars);
    vv=listVars{vvv};
    tmp1=read_nctiles([myenv.nctilesdir vv '/' vv],vv);
    tmp1=mean(tmp1,4);
    tmp1(mygrid.mskC==0)=NaN;
    eval([vv '=tmp1;']);
end;

UVELMASS=UVELMASS.*mygrid.mskW;
VVELMASS=VVELMASS.*mygrid.mskS;

if myenv.verbose>0; gcmfaces_msg('* call calc_barostream : comp. barotropic stream function');end;
[fldBAR]=calc_barostream(UVELMASS,VVELMASS);
if myenv.verbose>0; gcmfaces_msg('* call calc_overturn : comp. overturning stream function');end;
[gloOV]=calc_overturn(UVELMASS,VVELMASS);
if myenv.verbose>0; gcmfaces_msg('* call calc_transports : comp. transects transports');end;
[fldTRANSPORTS]=1e-6*calc_transports(UVELMASS,VVELMASS,mygrid.LINES_MASKS,{'dh','dz'});
if myenv.verbose>0; gcmfaces_msg('* call calc_MeridionalTransport : comp. meridional seawater transport');end;
[gloMT_FW]=1e-6*calc_MeridionalTransport(UVELMASS,VVELMASS,1);

if myenv.verbose>0; gcmfaces_msg('* load tracer and transports fields');end;
listVars={'THETA','SALT','ADVx_TH','ADVy_TH','ADVx_SLT','ADVy_SLT'};
listVars={listVars{:},'DFxE_TH','DFyE_TH','DFxE_SLT','DFyE_SLT'};
listDiags={listDiags{:},'fldTzonmean','fldSzonmean','gloMT_H','gloMT_SLT'};
for vvv=1:length(listVars);
    vv=listVars{vvv};
    tmp1=read_nctiles([myenv.nctilesdir vv '/' vv],vv);
    tmp1=mean(tmp1,4);
    tmp1(mygrid.mskC==0)=NaN;
    eval([vv '=tmp1;']);
end;

if myenv.verbose>0; gcmfaces_msg('* call calc_zonmean_T : comp. zonal mean temperature');end;
[fldTzonmean]=calc_zonmean_T(THETA);
if myenv.verbose>0; gcmfaces_msg('* call calc_zonmean_T : comp. zonal mean salinity');end;
[fldSzonmean]=calc_zonmean_T(SALT);

if myenv.verbose>0; gcmfaces_msg('* call calc_MeridionalTransport : comp. meridional heat transport');end;
tmpU=(ADVx_TH+DFxE_TH); tmpV=(ADVy_TH+DFyE_TH);
[gloMT_H]=1e-15*4e6*calc_MeridionalTransport(tmpU,tmpV,0);
if myenv.verbose>0; gcmfaces_msg('* call calc_MeridionalTransport : comp. meridional salt transport');end;
tmpU=(ADVx_SLT+DFxE_SLT); tmpV=(ADVy_SLT+DFyE_SLT);
[gloMT_SLT]=1e-6*calc_MeridionalTransport(tmpU,tmpV,0);

%format output:
for ddd=1:length(listDiags);
    dd=listDiags{ddd};
    eval(['diags.' dd '=' dd ';']);
end;

if myenv.verbose>0;
    gcmfaces_msg('*** leaving example_transports');
    gcmfaces_msg('===============================================','');
end;

