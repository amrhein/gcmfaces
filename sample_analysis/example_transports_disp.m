function []=example_transports_disp(diags);
%object:    display the result of example_transports
%inputs:    diags (the result of example_transports)

gcmfaces_global;
if myenv.verbose>0;
    gcmfaces_msg('===============================================');
    gcmfaces_msg(['*** entering example_transports_disp ' ...
        'that will display (as a figure or as text) the results of example_transports'],'');
end;


%%%%%%%%%%%%%%%
%display diags:
%%%%%%%%%%%%%%%

%barotropic streamfunction:
[X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,diags.fldBAR);
cc=[[-80:10:30] [40:40:200]];
figureL; set(gca,'FontSize',14);
pcolor(X,Y,FLD); axis([-180 180 -90 90]); 
set(gcf,'Renderer','zbuffer'); shading interp; 
gcmfaces_cmap_cbar(cc); title('Horizontal Stream Function');

%meridional streamfunction:
X=mygrid.LATS*ones(1,length(mygrid.RF)); Y=ones(length(mygrid.LATS),1)*mygrid.RF';
FLD=diags.gloOV; FLD(FLD==0)=NaN;
cc=[[-50:10:-30] [-24:3:24] [30:10:50]];
figureL; set(gca,'FontSize',14);
pcolor(X,Y,FLD); axis([-90 90 -6000 0]); 
set(gcf,'Renderer','zbuffer'); shading interp;
gcmfaces_cmap_cbar(cc); title('Meridional Stream Function');

if myenv.verbose>0; gcmfaces_msg('* call disp_transport : print and/or plot transports');end;

%Bering Strait and Arctic/Atlantic exchanges:
if length(diags.listTimes)>1; figureL; end;
transpList=[1 8:12];
rangeList=[[-1 3];[-3 1];[-6 2];[-3 9];[-9 3];[-0.5 0.5]];
for iii=1:length(transpList);
    if length(diags.listTimes)>1; subplot(3,2,iii); end;
    ylim=rangeList(iii,:);
    %
    ii=transpList(iii);
    trsp=diags.fldTRANSPORTS(ii,:)';
    txt=[mygrid.LINES_MASKS(ii).name ' (>0 to Arctic)'];
    disp_transport(trsp,diags.listTimes,txt,{'ylim',ylim});
end;

%Drake, ACC etc:
if length(diags.listTimes)>1; figureL; end;
transpList=[13 20 19 18];
rangeList=[[120 200];[120 200];[-40 10];[120 200]];
for iii=1:length(transpList);
    if length(diags.listTimes)>1; subplot(3,2,iii); end;
    ylim=rangeList(iii,:);
    %
    ii=transpList(iii);
    trsp=diags.fldTRANSPORTS(ii,:)';
    txt=[mygrid.LINES_MASKS(ii).name ' (>0 to the West)'];
    disp_transport(trsp,diags.listTimes,txt,{'ylim',ylim});
end;

%%%%

X=mygrid.LATS*ones(1,length(mygrid.RC)); Y=ones(length(mygrid.LATS),1)*mygrid.RC';
FLD=diags.fldTzonmean; FLD(FLD==0)=NaN;
cc=[-3:2:30];
figureL; set(gca,'FontSize',14);
pcolor(X,Y,FLD); axis([-90 90 -6000 0]); 
set(gcf,'Renderer','zbuffer'); shading interp;
gcmfaces_cmap_cbar(cc); title('zonal mean THETA');

X=mygrid.LATS*ones(1,length(mygrid.RC)); Y=ones(length(mygrid.LATS),1)*mygrid.RC';
FLD=diags.fldSzonmean; FLD(FLD==0)=NaN;
cc=[32:0.2:36];
figureL; set(gca,'FontSize',14);
pcolor(X,Y,FLD); axis([-90 90 -6000 0]); 
set(gcf,'Renderer','zbuffer'); shading interp;
gcmfaces_cmap_cbar(cc); title('zonal mean SALT');

FLD=diags.gloMT_H; FLD(FLD==0)=NaN;
figureL; set(gca,'FontSize',14);
plot(mygrid.LATS,FLD); 
title('Meridional Heat Transport (in PW)');

FLD=diags.gloMT_FW; FLD(FLD==0)=NaN;
figureL; set(gca,'FontSize',14);
plot(mygrid.LATS,FLD); 
title('Meridional Sea Water Transport (in Sv)');

FLD=diags.gloMT_SLT; FLD(FLD==0)=NaN;
figureL; set(gca,'FontSize',14);
plot(mygrid.LATS,FLD); 
title('Meridional Salt Transport (in psu Sv)');

%%%%

if myenv.verbose>0;
    gcmfaces_msg('*** leaving example_transports_disp');
    gcmfaces_msg('===============================================','');
end;
