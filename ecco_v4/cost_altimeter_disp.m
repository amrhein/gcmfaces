function []=cost_altimeter_disp(dirMat,choicePlot,suf,dirTex,nameTex);
%object:	plot the various sea level statistics
%		(std model-obs, model, obs, leading to cost function terms)
%inputs:	dirMat is the model run directory
%		choicePlot is 1 (rms) 2 (prior uncertainty) or 3 (cost)
%               suf is 'modMobs', 'obs' or 'mod'
%optional:      dirTex is the directory where tex and figures files are created
%                 (if not specified then display all results to screen instead)
%               nameTex is the tex file name (default : 'myPlots')

gcmfaces_global;

%backward compatibility test
test1=~isempty(dir([dirMat 'basic_diags_ecco_mygrid.mat']));
test2=~isempty(dir([dirMat 'diags_grid_parms.mat']));
if ~test1&~test2;
    error('missing diags_grid_parms.mat')
elseif test2;
    nameGrid='diags_grid_parms.mat';
else;
    nameGrid='basic_diags_ecco_mygrid.mat';
end;

%here we always reload the grid from dirMat to make sure the same one is used throughout
eval(['load ' dirMat nameGrid ';']);

%determine if and where to create tex and figures files
dirMat=[dirMat '/'];
if isempty(who('dirTex'));
    addToTex=0;
else;
    if ~ischar(dirTex); error('mis-specified dirTex'); end;
    addToTex=1;
    if isempty(who('nameTex')); nameTex='myPlots'; end;
    fileTex=[dirTex nameTex '.tex'];
end;

%%%%%%%%%%%%%%%
%define pathes:
%%%%%%%%%%%%%%%

if isempty(dirMat); dirMat=[dirModel 'mat/']; else; dirMat=[dirMat '/']; end;
runName=pwd; tmp1=strfind(runName,filesep); runName=runName(tmp1(end)+1:end);

%%%%%%%%%%%%%%%%%
%do computations:
%%%%%%%%%%%%%%%%%

if isdir([dirMat 'cost/']); dirMat=[dirMat 'cost/']; end;

eval(['load ' dirMat 'cost_altimeter_' suf '.mat myflds;']);

%mask for plotting
if isfield(myflds,'msk_100pts');
    MSK=myflds.msk_100pts;
else;
    MSK=mygrid.mskC(:,:,1);
end;
%hack to omit any masking : MSK(:)=1;

if strcmp(suf,'modMobs'); tit='modeled-observed';
elseif strcmp(suf,'obs'); tit='observed';
elseif strcmp(suf,'mod'); tit='modeled';
else; error('unknown field');
end

if choicePlot==1;  	tit=[tit ' log(variance)']; uni='(m$^2$)';
elseif choicePlot==2;	tit=' log(prior error variance)'; uni='(m$^2$)';
else; 			tit=[tit ' cost']; uni='';
end;

if choicePlot==1;%rms
    
    if strcmp(suf,'modMobs');
        cc=[-0.4:0.05:-0.25 -0.2:0.03:-0.05 -0.03:0.01:0.03 0.05:0.03:0.2 0.25:0.05:0.4];
        figure; m_map_gcmfaces(100*MSK.*myflds.dif_mdt,0,{'myCaxis',100*cc}); drawnow;
        myCaption={'mean dynamic topography misfit (cm)'};
        if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
    end;
    
    cc=[-4:0.2:-1]-0.5;
    figure; m_map_gcmfaces(log10(MSK.*(myflds.rms_sladiff_smooth.^2)),0,{'myCaxis',cc}); drawnow;
    myCaption={tit,'-- sea level anomaly ',uni,' -- large space/time scales'};
    if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
    
    cc=[-4:0.2:-1];
    figure; m_map_gcmfaces(log10(MSK.*(myflds.rms_sladiff_point.^2)),0,{'myCaxis',cc}); drawnow;
    myCaption={tit,'-- sea level anomaly ',uni,' -- pointwise'};
    if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
    
    if 0;
        figure; m_map_gcmfaces(log10(MSK.*(myflds.rms_sladiff_point35d.^2)),0,{'myCaxis',cc}); drawnow;
        myCaption={tit,'-- sea level anomaly ',uni,' -- large time scales'};
        if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
        
        figure; m_map_gcmfaces(log10(MSK.*(myflds.rms_sladiff_35dMsmooth.^2)),0,{'myCaxis',cc}); drawnow;
        myCaption={tit,'-- sea level anomaly ',uni,' -- point35d minus lsc'};
        if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
        
        figure; m_map_gcmfaces(log10(MSK.*(myflds.rms_sladiff_pointMpoint35d.^2)),0,{'myCaxis',cc}); drawnow;
        myCaption={tit,'-- sea level anomaly ',uni,' -- point minus point35d'};
        if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
    end;
    
elseif choicePlot==2;%uncertainty fields
    
    cc=[0:0.005:0.02 0.03:0.01:0.05 0.06:0.02:0.1 0.14:0.03:0.2 0.25:0.05:0.4];
    figure; m_map_gcmfaces(100*myflds.sig_mdt,0,{'myCaxis',100*cc}); drawnow;
    myCaption={'mean dynamic topography prior uncertainty (cm)'};
    if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
    
    cc=[-4:0.2:-1]-0.5;
    figure; m_map_gcmfaces(log10(MSK.*(myflds.sig_sladiff_smooth.^2)),0,{'myCaxis',cc}); drawnow;
    myCaption={tit,'-- sea level anomaly ',uni,' -- large space/time scales'};
    if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
    
    cc=[-4:0.2:-1];
    figure; m_map_gcmfaces(log10(MSK.*(myflds.sig_sladiff_point.^2)),0,{'myCaxis',cc}); drawnow;
    myCaption={tit,'-- sea level anomaly ',uni,' -- pointwise'};
    if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
    
else;%cost
    
    cc=[0:0.005:0.02 0.03:0.01:0.05 0.06:0.02:0.1 0.14:0.03:0.2 0.25:0.05:0.4]*100;
    
    figure; m_map_gcmfaces(MSK.*((myflds.dif_mdt.^2)./(myflds.sig_mdt.^2)),0,{'myCaxis',cc}); drawnow;
    myCaption={tit,'-- mean dynamic topography ',uni};
    if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
    
    figure; m_map_gcmfaces(MSK.*((myflds.rms_sladiff_smooth.^2)./(myflds.sig_sladiff_smooth.^2)),0,{'myCaxis',cc}); drawnow;
    myCaption={tit,'-- sea level anomaly ',uni,' -- large space/time scales'};
    if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
    
    figure; m_map_gcmfaces(MSK.*((myflds.rms_sladiff_point.^2)./(myflds.sig_sladiff_point.^2)),0,{'myCaxis',cc}); drawnow;
    myCaption={tit,'-- sea level anomaly ',uni,' -- pointwise'};
    if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
    
end;





