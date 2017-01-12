function []=plot_one_field(choiceV3orV4,plotMethod);
%object:    illustrate several methods to plot fields within gcmfaces.
%inputs:    choiceV3orV4 ('v3' or 'v4') selects the sample GRID
%           plotMethod selects the plotting method
%               0) all, one after the other
%               1) plot a gcmfaces object face by face
%               2) use convert2array then plot in array coordinates.
%               3) use convert2pcol then plot in lat-lon coordinates.
%               4) use the gcmfaces front end to m_map.

input_list_check('plot_one_field',nargin);

gcmfaces_global;
if myenv.verbose>0;
    gcmfaces_msg('===============================================');
    gcmfaces_msg(['*** entering plot_one_field that displays ' ... 
        'a gridded field of ocean depth, first ' ... 
        'in gcmfaces and array formats. Then in geographic coordinates ' ...
        'using pcolor and in various projections using m_map (if avail. in path). '],'');
end;

%%%%%%%%%%%%%%%%%
%load parameters:
%%%%%%%%%%%%%%%%%

gcmfaces_global;
myenv.nctiles=1;
if strcmp(choiceV3orV4,'v4');
    nF=5;
    fileFormat='compact';
    dir0=fullfile([myenv.gcmfaces_dir '../']);
    dirGrid=fullfile(dir0,'GRID/');
else;
    error('files are missing');
    nF=1;
    dir0=fullfile(myenv.gcmfaces_dir,'sample_input/');
    fileFormat='straight';
    dirGrid=fullfile(dir0,'grid_occa/');
end;

mygrid=[]; grid_load(dirGrid,nF,fileFormat,1);

%%%%%%%%%%%
%get field:
%%%%%%%%%%%

fld=mygrid.Depth; fld(fld==0)=NaN;
cc=[[0:0.05:0.5] [0.6 0.75 1 1.25]]*1e4; myCmap='gray';

%%%%%%%%%%%%
%plot field:
%%%%%%%%%%%%

if plotMethod==0|plotMethod==1;
    if myenv.verbose>0; gcmfaces_msg('* gcmfaces format display -- face by face.'); end;
    if nF==1;
        figureL; imagescnan(fld{1}','nancolor',[1 1 1]*0.8); axis xy; cb=gcmfaces_cmap_cbar(cc,{'myCmap',myCmap}); delete(cb);
    elseif nF==5;
        figureL;
        subplot(3,3,7); imagescnan(fld{1}','nancolor',[1 1 1]*0.8); axis xy; cb=gcmfaces_cmap_cbar(cc,{'myCmap',myCmap}); delete(cb);
        tmp1=axis; tmp2=text(tmp1(2)/2,tmp1(4)/2,'1','FontSize',32,'Color','r','Rotation',0);
        subplot(3,3,8); imagescnan(fld{2}','nancolor',[1 1 1]*0.8); axis xy; cb=gcmfaces_cmap_cbar(cc,{'myCmap',myCmap}); delete(cb);
        tmp1=axis; tmp2=text(tmp1(2)/2,tmp1(4)/2,'2','FontSize',32,'Color','r','Rotation',0);
        subplot(3,3,5); imagescnan(fld{3}','nancolor',[1 1 1]*0.8); axis xy; cb=gcmfaces_cmap_cbar(cc,{'myCmap',myCmap}); delete(cb);
        tmp1=axis; tmp2=text(tmp1(2)/2,tmp1(4)/2,'3','FontSize',32,'Color','r','Rotation',0);
        subplot(3,3,6); imagescnan(fld{4}','nancolor',[1 1 1]*0.8); axis xy; cb=gcmfaces_cmap_cbar(cc,{'myCmap',myCmap}); delete(cb);
        tmp1=axis; tmp2=text(tmp1(2)/2,tmp1(4)/2,'4','FontSize',32,'Color','r','Rotation',0);
        subplot(3,3,3); imagescnan(fld{5}','nancolor',[1 1 1]*0.8); axis xy; cb=gcmfaces_cmap_cbar(cc,{'myCmap',myCmap}); delete(cb);
        tmp1=axis; tmp2=text(tmp1(2)/2,tmp1(4)/2,'5','FontSize',32,'Color','r','Rotation',0);
    elseif nF==6;
        error('face by face plot for cude : not yet implemented');
    end;
end;

if plotMethod==0|plotMethod==2;
    if myenv.verbose>0; gcmfaces_msg('* array format display -- all faces concatenated somehow.'); end;
    figureL;
    FLD=convert2array(fld);
    imagescnan(FLD','nancolor',[1 1 1]*0.8); axis xy; cb=gcmfaces_cmap_cbar(cc,{'myCmap',myCmap}); %delete(cb); 
    for ff=1:mygrid.nFaces; 
        tmp0=0*mygrid.XC;
        tmp1=round(size(tmp0{ff})/2); 
        tmp0{ff}(tmp1(1),tmp1(2))=1; 
        tmp0=convert2array(tmp0); [ii,jj]=find(tmp0==1);
        if ff==3; ang=90; elseif ff>3; ang=-90; else; ang=0; end;
        hold on; text(ii,jj,num2str(ff),'FontSize',32,'Color','r','Rotation',ang);
    end;
end;

if plotMethod==0|plotMethod==3;
    if myenv.verbose>0; gcmfaces_msg('* geographical display -- using pcolor.'); end;
    figureL;
    [X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,fld); pcolor(X,Y,FLD);
    if ~isempty(find(X>359)); axis([0 360 -90 90]); else; axis([-180 180 -90 90]); end;
    shading flat; cb=gcmfaces_cmap_cbar(cc,{'myCmap',myCmap}); delete(cb);  xlabel('longitude'); ylabel('latitude');
end;

if plotMethod==0|plotMethod==4;
    if myenv.verbose>0; gcmfaces_msg('* geographical display -- using m_map.'); end;
    if ~isempty(which('m_proj'));
        figureL; m_map_gcmfaces(fld,0,{'myCaxis',cc},{'myCmap',myCmap});
    elseif myenv.verbose;
        fprintf('if you want to use m_map_gcmfaces, you will need to have m_map in your path\n');
    end;
end;

if plotMethod==0|plotMethod==5;
    if myenv.verbose>0; gcmfaces_msg('* geographical display -- on a sphere.'); end;
    gcmfaces_sphere(fld,cc,[],'N',3);
end;

%test case for depthStretchPlot:
if plotMethod==0|plotMethod==6;
    if myenv.verbose>0; gcmfaces_msg('* section display -- using strecthed vertical coord.'); end;
    x=ones(length(mygrid.RC),1)*[1:200]; z=mygrid.RC*ones(1,200); c=sin(z/2000*pi).*cos(x/50*pi);
    figureL; 
    subplot(1,2,1); pcolor(x,z,c); shading flat; title('standard depth display');
    subplot(1,2,2); depthStretchPlot('pcolor',{x,z,c}); shading flat; title('stretched depth display');
%     subplot(1,2,1); plot(c(:,1),z(:,1)); title('standard depth display');
%     subplot(1,2,2); depthStretchPlot('plot',{c(:,1),z(:,1)}); title('stretched depth display');
end;

if myenv.verbose>0;
    gcmfaces_msg('*** leaving plot_one_field');
    gcmfaces_msg('===============================================');
end;

