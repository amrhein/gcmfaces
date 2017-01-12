function []=example_budget(choiceV3orV4);
%object:    illustrates budget analysis
%inputs:    choiceV3orV4 ('v3' or 'v4') selects the sample GRID

gcmfaces_global;
if myenv.verbose>0;
    gcmfaces_msg('===============================================');
    gcmfaces_msg(['*** entering example_budget ' ...
        'that will load a set of budget terms from file ' ...
        'and compute a series of derived diagnostics '],'');
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
    myenv.nctilesdir=fullfile(dir0,'release1/nctiles_budget/');
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

%%%%%%%%%%%%%%%%%%
%main computation:
%%%%%%%%%%%%%%%%%%

%select budget of interest
nameBudg='budgMo';

%load budget terms
fileName=[myenv.nctilesdir nameBudg filesep nameBudg(1:6)];
tend=read_nctiles(fileName,'tend');
trU=read_nctiles(fileName,'trU');
trV=read_nctiles(fileName,'trV');
trWtop=read_nctiles(fileName,'trWtop');
trWbot=read_nctiles(fileName,'trWbot');

%load dt (time increments) vector (not a gcmfaces object)
ncload([fileName '.0001.nc'],'dt');

%get budget descitption and units
ncid = netcdf.open([fileName '.0001.nc'],'NC_NOWRITE');
descr = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'description');
varid = netcdf.inqVarID(ncid,'tend');
units = netcdf.getAtt(ncid,varid,'units');
netcdf.close(ncid);

%define northern hemisphere as domain of integration
nameMask='Northern Hemisphere';
mask=mygrid.mskC(:,:,1).*(mygrid.YC>0); 
areaMask=mygrid.RAC.*mask;

%edit plot title accordingly
tmp1=strfind(descr,'-- ECCO v4');
descr=[descr(1:tmp1-1) 'for: ' nameMask];

%compute northern hemisphere integrals
budg.tend=NaN*dt;
budg.hconv=NaN*dt;
budg.zconv=NaN*dt;
for tt=1:length(dt);
    %compute flux convergence
    hconv=calc_UV_conv(trU(:,:,tt),trV(:,:,tt));
    zconv=trWtop(:,:,tt)-trWbot(:,:,tt);
    %compute sum over domain
    budg.tend(tt)=nansum(tend(:,:,tt).*mask)/nansum(areaMask);
    budg.hconv(tt)=nansum(hconv.*mask)/nansum(areaMask);
    budg.zconv(tt)=nansum(zconv.*mask)/nansum(areaMask);
end;

%display result
figureL;
subplot(3,1,1); set(gca,'FontSize',12);
plot(cumsum(dt.*budg.tend));
grid on; xlabel('month'); ylabel([units '.s']);
legend('content anomaly'); title(descr);
subplot(3,1,2); set(gca,'FontSize',12);
plot(cumsum(dt.*budg.tend)); hold on;
plot(cumsum(dt.*budg.hconv),'r'); 
plot(cumsum(dt.*budg.zconv),'g');
grid on; xlabel('month'); ylabel([units '.s']);
legend('content anomaly','horizontal convergence','vertical convergence');
subplot(3,1,3); set(gca,'FontSize',12);
plot(cumsum(dt.*(budg.tend-budg.hconv-budg.zconv))); 
grid on; xlabel('month'); ylabel([units '.s']);
legend('budget residual');

