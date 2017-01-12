function []=cost_sst(dirModel,dirMat,doComp,dirTex,nameTex);
%object:        compute cost function term for sst data
%inputs:        dimodel is the model directory
%               dirMat is the directory where diagnozed .mat files will be saved
%                     -> set it to '' to use the default [dirModel 'mat/']
%               doComp is a switch (1->compute; 0->display)
%optional:      dirTex is the directory where tex and figures files are created
%                 (if not specified then display all results to screen instead)
%               nameTex is the tex file name (default : 'myPlots')

if isempty(dirMat); dirMat=[dirModel 'mat/']; else; dirMat=[dirMat '/']; end;
if isempty(dir(dirMat));     mkdir([dirMat]); end;

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

if doComp;

%grid, params and inputs

gcmfaces_global; global myparms;
if ~isfield(mygrid,'XC'); grid_load('./GRID/',5,'compact'); end;
if ~isfield(mygrid,'LATS_MASKS'); gcmfaces_lines_zonal; end;
if isfield(myparms,'yearFirst'); yearFirst=myparms.yearFirst; yearLast=myparms.yearLast;
else; yearFirst=1992; yearLast=2011;
end;

dirData='/net/nares/raid11/ecco-shared/ecco-version-4/input/';
subdirReynolds='input_sst_reynolds_v4/'; nameReynolds='reynolds_oiv2_r1';
subdirRemss='input_sst_remss_v4/'; nameRemss='tmi_amsre_oisst_r1';
subdirL2p='input_sst_amsre_daily/'; nameL2p='tmi_amsre_l2p_r1';

if 0;%old versions
nameReynolds='v4_reynolds_monthly';
nameRemss='v4_TMI_AMSRE_sst_monav';
nameL2p='amsre_r2';
yearLast=2010;
end;

dirErr='/net/nares/raid11/ecco-shared/ecco-version-4/input/input_all_from_pleiades/';
fld_err=read2memory([dirErr '/sigma_half.bin'],[90 1170]);
fld_err=convert2gcmfaces(fld_err);
fld_w=fld_err.^-2;

fileModel=dir([dirModel 'barfiles/tbar*data']); fileModel=['barfiles/' fileModel.name];

%computational loop
mod_m_rey=repmat(NaN*mygrid.mskC(:,:,1),[1 1 12*(yearLast-yearFirst+1)]);
mod_m_remss=repmat(NaN*mygrid.mskC(:,:,1),[1 1 12*(yearLast-yearFirst+1)]);
%
zm_mod=repmat(ones(length(mygrid.LATS),1),[1 12*(yearLast-yearFirst+1)]);
zm_rey=repmat(ones(length(mygrid.LATS),1),[1 12*(yearLast-yearFirst+1)]);
zm_mod_m_rey=repmat(ones(length(mygrid.LATS),1),[1 12*(yearLast-yearFirst+1)]);
zm_remss=repmat(ones(length(mygrid.LATS),1),[1 12*(yearLast-yearFirst+1)]);
zm_mod_m_remss=repmat(ones(length(mygrid.LATS),1),[1 12*(yearLast-yearFirst+1)]);
%
for ycur=yearFirst:yearLast;
fprintf(['starting ' num2str(ycur) '\n']);
tic;
for mcur=1:12;
  %load Reynolds SST
  file0=[dirData subdirReynolds nameReynolds '_' num2str(ycur)];
  if ~isempty(dir(file0)); fld_rey=read_bin(file0,mcur,0); else; fld_rey=NaN*mygrid.mskC(:,:,1); end;
  fld_rey(find(fld_rey==0))=NaN;
  fld_rey(find(fld_rey<-99))=NaN;

  %load Remss SST
  file0=[dirData subdirRemss nameRemss '_' num2str(ycur)];
  if ~isempty(dir(file0)); fld_remss=read_bin(file0,mcur,0); else; fld_remss=NaN*mygrid.mskC(:,:,1); end;
  fld_remss(find(fld_remss==0))=NaN;
  fld_remss(find(fld_remss<-99))=NaN;

  %load model SST
  mm=(ycur-yearFirst)*12+mcur;
  fld_mod=read_bin([dirModel fileModel],mm,1).*mygrid.mskC(:,:,1);

  %store misfit maps 
  mod_m_rey(:,:,mm)=fld_mod-fld_rey;
  mod_m_remss(:,:,mm)=fld_mod-fld_remss;

  %comp zonal means
  zm_mod(:,mm)=calc_zonmean_T(fld_mod.*mygrid.mskC(:,:,1));

  msk=mygrid.mskC(:,:,1); msk(fld_rey==0)=NaN;
  zm_rey(:,mm)=calc_zonmean_T(fld_rey.*msk);
  zm_mod_m_rey(:,mm)=calc_zonmean_T(fld_mod-fld_rey.*msk);

  msk=mygrid.mskC(:,:,1); msk(fld_remss==0)=NaN;
  zm_remss(:,mm)=calc_zonmean_T(fld_remss.*msk);
  zm_mod_m_remss(:,mm)=calc_zonmean_T(fld_mod-fld_remss.*msk);
end;
toc;
end;

%compute rms maps
rms_to_rey=sqrt(nanmean(mod_m_rey.^2,3));
rms_to_remss=sqrt(nanmean(mod_m_remss.^2,3));

if ~isdir([dirMat 'cost/']); mkdir([dirMat 'cost/']); end;
eval(['save ' dirMat '/cost/cost_sst.mat fld_err rms_* zm_*;']);

else;%display previously computed results

global mygrid;

if isdir([dirMat 'cost/']); dirMat=[dirMat 'cost/']; end;

eval(['load ' dirMat '/cost_sst.mat;']);

if ~isempty(who('fld_rms')); 
  figure; m_map_gcmfaces(fld_rms,0,{'myCaxis',[0:0.2:1.2 1.5:0.5:3 4:1:6 8 10]/2});
  myCaption={'modeled-observed rms -- sea surface temperature (K)'};
  if addToTex; write2tex(fileTex,2,myCaption,gcf); end;
else;

  figure; m_map_gcmfaces(rms_to_rey,0,{'myCaxis',[0:0.2:1.2 1.5:0.5:3 4:1:6 8 10]/2});
  myCaption={'modeled-Reynolds rms -- sea surface temperature (K)'};
  if addToTex; write2tex(fileTex,2,myCaption,gcf); end;

  figure; m_map_gcmfaces(rms_to_remss,0,{'myCaxis',[0:0.2:1.2 1.5:0.5:3 4:1:6 8 10]/2});
  myCaption={'modeled-REMSS rms -- sea surface temperature (K)'};
  if addToTex; write2tex(fileTex,2,myCaption,gcf); end;

  ny=size(zm_rey,2)/12;
  [xx,yy]=meshgrid(1992+([1:ny*12]-0.5)/12,mygrid.LATS);

  figureL;
  obs=zm_rey; obsCycle=sst_cycle(obs); 
  mod=zm_rey+zm_mod_m_rey; modCycle=sst_cycle(mod);
  mis=zm_mod_m_rey; misCycle=sst_cycle(mis);
  subplot(3,1,1); pcolor(xx,yy,obs-obsCycle); shading flat; caxis([-1 1]*1); colorbar;
  set(gca,'FontSize',14); set(gca,'XTick',[]); ylabel('latitude'); 
  title('Reynolds sst anomaly');
  subplot(3,1,2); pcolor(xx,yy,mod-modCycle); shading flat; caxis([-1 1]*1); colorbar;
  set(gca,'FontSize',14); set(gca,'XTick',[]); ylabel('latitude');
  title('ECCO sst anomaly');
  subplot(3,1,3); pcolor(xx,yy,mis-misCycle); shading flat; caxis([-1 1]*1); colorbar;
  set(gca,'FontSize',14); ylabel('latitude'); title('ECCO-Reynolds sst misfit');
  myCaption={'ECCO and Reynolds zonal mean sst anomalies (K)'};
  if addToTex; write2tex(fileTex,2,myCaption,gcf); end;

  figureL;
  obs=zm_remss; obsCycle=sst_cycle(obs);
  mod=zm_remss+zm_mod_m_remss; modCycle=sst_cycle(mod);
  mis=zm_mod_m_remss; misCycle=sst_cycle(mis);
  subplot(3,1,1); pcolor(xx,yy,obs-obsCycle); shading flat; caxis([-1 1]*1); colorbar;
  set(gca,'FontSize',14); set(gca,'XTick',[]); ylabel('latitude');
  title('REMSS sst anomaly');
  subplot(3,1,2); pcolor(xx,yy,mod-modCycle); shading flat; caxis([-1 1]*1); colorbar;
  set(gca,'FontSize',14); set(gca,'XTick',[]); ylabel('latitude');
  title('ECCO sst anomaly');
  subplot(3,1,3); pcolor(xx,yy,mis-misCycle); shading flat; caxis([-1 1]*1); colorbar;
  set(gca,'FontSize',14); ylabel('latitude'); title('ECCO-REMSS sst misfit');
  myCaption={'ECCO and REMSS zonal mean sst anomalies (K)'};
  if addToTex; write2tex(fileTex,2,myCaption,gcf); end;

end;

end;

function [zmCycle]=sst_cycle(zmIn);

ny=size(zmIn,2)/12;
zmCycle=NaN*zeros(179,12);
for mm=1:12;
zmCycle(:,mm)=nanmean(zmIn(:,mm:12:ny*12),2);
end;
zmCycle=repmat(zmCycle,[1 ny]);

