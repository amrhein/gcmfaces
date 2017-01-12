
%select kBudget:
if ~isempty(setDiagsParams);
  kBudget=setDiagsParams{1};
else;
  kBudget=1;
end;

%override default file name:
%---------------------------
tmp1=setDiags;
if kBudget>1;
    tmp1=sprintf('E%02i',kBudget);
end;
fileMat=['diags_set_' tmp1];

if userStep==1;%diags to be computed
    listDiags=['zm_area zm_vol_ocn zm_vol_tot zm_vol_ice'];
    listDiags=[listDiags ' zm_heat_ocn zm_heat_tot zm_heat_ice zm_heat_ocn_diff'];
    listDiags=[listDiags ' zm_salt_ocn zm_salt_tot zm_salt_ice zm_salt_ocn_diff'];
elseif userStep==2;%input files and variables
    listFlds={    'ETAN','SIheff','SIhsnow','THETA   ','SALT    ','PHIBOT'};
    listFlds={listFlds{:},'SIatmFW ','oceFWflx','SItflux','TFLUX','SFLUX','oceSPflx','SRELAX'};
    listFlds={listFlds{:},'oceQnet ','SIatmQnt','SIaaflux','SIsnPrcp','SIacSubl'};
    listFlds={listFlds{:},'TRELAX','WTHMASS','WSLTMASS','oceSflux','oceQsw','oceSPtnd'};
    if kBudget>1;
        listFlds={listFlds{:},'ADVr_TH','DFrE_TH','DFrI_TH','ADVr_SLT','DFrE_SLT','DFrI_SLT','WVELMASS'};
    end;
    listFlds={listFlds{:},'UVELMASS','VVELMASS','AB_gT','AB_gS'};
    listFlds={listFlds{:},'ADVx_TH ','ADVy_TH ','DFxE_TH ','DFyE_TH '};
    listFlds={listFlds{:},'ADVx_SLT','ADVy_SLT','DFxE_SLT','DFyE_SLT'};
    listFlds={listFlds{:},'ADVxHEFF','ADVyHEFF','DFxEHEFF','DFyEHEFF'};
    listFlds={listFlds{:},'ADVxSNOW','ADVySNOW','DFxESNOW','DFyESNOW'};
    listFldsNames=deblank(listFlds);
    %
    listFiles={'rate_budg2d_snap_set1','budg2d_hflux_set1','budg2d_zflux_set1','budg2d_zflux_set2'};
    if kBudget==1;
        listFiles={listFiles{:},'rate_budg2d_snap_set2','budg2d_hflux_set2'};
    else;
        tmp1=sprintf('rate_budg2d_snap_set3_%02i',kBudget);
        tmp2=sprintf('budg2d_zflux_set3_%02i',kBudget);
        tmp3=sprintf('budg2d_hflux_set3_%02i',kBudget);
        listFiles={listFiles{:},tmp1,tmp2,tmp3};
    end;
    listSubdirs={[dirMat 'BUDG/' ],[dirMat '../BUDG/' ],[dirModel 'diags/BUDG/'],[dirModel 'diags/']};
elseif userStep==3;%computational part;

    %preliminary tests
    test1=isempty(dir([dirModel 'diags/BUDG/budg2d_snap_set1*']));
    test2=isempty(dir([dirMat 'BUDG/rate_budg2d_snap_set1*']))&...
          isempty(dir([dirMat '../BUDG/rate_budg2d_snap_set1*']));

    if (strcmp(setDiags,'E')&test1&test2);
        fprintf('\n abort : global and regional budgets, due to missing \n');
        fprintf(['\n   ' dirModel 'diags/BUDG/budg2d_snap_set1* \n']);
        return;
    end;

    if (strcmp(setDiags,'E')&test2);
        fprintf('\n abort : global and regional budgets, due to missing \n');
        fprintf(['\n   ' dirModel 'diags/BUDG/rate_budg2d_snap_set1* \n']);
        return;
    end;
    
    %override default file name:
    %---------------------------
    tmp1=setDiags;
    if kBudget>1;
        tmp1=sprintf('E%02i',kBudget);
    end;
    fileMat=['diags_set_' tmp1 '_' num2str(tt) '.mat'];
    
    %fill in optional fields:
    %------------------------
    if isempty(who('TRELAX')); TRELAX=0*mygrid.XC; end;
    if isempty(who('SRELAX')); SRELAX=0*mygrid.XC; end;
    if isempty(who('AB_gT')); AB_gT=0*mygrid.XC; end;
    if isempty(who('AB_gS')); AB_gS=0*mygrid.XC; end;
    if isempty(who('oceSPtnd')); oceSPtnd=0*mygrid.XC; end;
    if isempty(who('oceSPflx')); oceSPflx=0*mygrid.XC; end;
    if isempty(who('PHIBOT')); PHIBOT=0*mygrid.XC; end;
    
    %=======MASS=========
    
    %compute mapped budget:
    %----------------------
    
    %mass = myparms.rhoconst * sea level
    contOCN=ETAN*myparms.rhoconst;
    contICE=(SIheff*myparms.rhoi+SIhsnow*myparms.rhosn);
    %for deep ocean layer :
    if kBudget>1&myparms.useNLFS<2;
        contOCN=0;
    elseif kBudget>1;%rstar case
        tmp1=mk3D(mygrid.DRF,mygrid.hFacC).*mygrid.hFacC;
        tmp2=sum(tmp1(:,:,kBudget:length(mygrid.RC)),3)./mygrid.Depth;
        contOCN=tmp2.*ETAN*myparms.rhoconst;
    end;
    %
    contTOT=contOCN+contICE;
    %vertical divergence (air-sea fluxes or vertical advection)
    zdivOCN=oceFWflx;
    zdivICE=SIatmFW-oceFWflx;
    %in virtual salt flux we omit :
    if ~myparms.useRFWF; zdivOCN=0*zdivOCN; end;
    
    %for deep ocean layer :
    if kBudget>1; zdivOCN=-WVELMASS*myparms.rhoconst; end;
    %
    zdivTOT=zdivOCN+zdivICE;
    %horizontal divergence (advection and ice diffusion)
    hdivOCN=myparms.rhoconst*calc_UV_conv(UVELMASS,VVELMASS,{'dh'}); %for 2D, already vertically integrated, fields
    tmpU=(myparms.rhoi*DFxEHEFF+myparms.rhosn*DFxESNOW+myparms.rhoi*ADVxHEFF+myparms.rhosn*ADVxSNOW);
    tmpV=(myparms.rhoi*DFyEHEFF+myparms.rhosn*DFyESNOW+myparms.rhoi*ADVyHEFF+myparms.rhosn*ADVySNOW);
    hdivICE=calc_UV_conv(tmpU,tmpV); %no dh needed here
    hdivTOT=hdivOCN+hdivICE;
    %bottom pressure for comparison:
    bp=myparms.rhoconst/9.81*PHIBOT;
    
    [zm_vol_tot,zm_area]=calc_budget_mean_zonal(contTOT,zdivTOT,hdivTOT);
    zm_vol_ocn=calc_budget_mean_zonal(contOCN,zdivOCN,hdivOCN);
    zm_vol_ice=calc_budget_mean_zonal(contICE,zdivICE,hdivICE);
    
    %=======HEAT=======
    
    contOCN=myparms.rcp*THETA-myparms.rcp*AB_gT;
    contICE=-myparms.flami*(SIheff*myparms.rhoi+SIhsnow*myparms.rhosn);
    contTOT=contOCN+contICE;
    %vertical divergence (air-sea fluxes or vertical adv/dif)
    zdivOCN=TFLUX;
    zdivICE=-(SItflux+TFLUX-TRELAX);
    %in linear surface we omit :
    if ~myparms.useNLFS; zdivOCN=zdivOCN-myparms.rcp*WTHMASS; end;
    %in virtual salt flux we omit :
    if ~myparms.useRFWF|~myparms.useNLFS; zdivICE=zdivICE+SIaaflux; end;
    %working approach for real fresh water (?) and virtual salt flux
    if 0; zdivICE=-oceQnet-SIatmQnt-myparms.flami*(SIsnPrcp-SIacSubl); end;
    %for deep ocean layer :
    if kBudget>1;
        zdivOCN=-(ADVr_TH+DFrE_TH+DFrI_TH)./mygrid.RAC*myparms.rcp;
        dd=mygrid.RF(kBudget); msk=mygrid.mskC(:,:,kBudget);
        swfrac=0.62*exp(dd/0.6)+(1-0.62)*exp(dd/20);
        if dd<-200; swfrac=0; end;
        zdivOCN=zdivOCN+swfrac*oceQsw;%.*msk;
    end;
    %
    zdivTOT=zdivOCN+zdivICE;
    %horizontal divergence (advection and diffusion)
    tmpU=myparms.rcp*(ADVx_TH+DFxE_TH); tmpV=myparms.rcp*(ADVy_TH+DFyE_TH);
    hdivOCN=calc_UV_conv(tmpU,tmpV);
    tmpU=-myparms.flami*(myparms.rhoi*DFxEHEFF+myparms.rhosn*DFxESNOW+myparms.rhoi*ADVxHEFF+myparms.rhosn*ADVxSNOW);
    tmpV=-myparms.flami*(myparms.rhoi*DFyEHEFF+myparms.rhosn*DFyESNOW+myparms.rhoi*ADVyHEFF+myparms.rhosn*ADVySNOW);
    hdivICE=calc_UV_conv(tmpU,tmpV); %no dh needed here
    hdivTOT=hdivOCN+hdivICE;
    
    zm_heat_tot=calc_budget_mean_zonal(contTOT,zdivTOT,hdivTOT);
    zm_heat_ocn=calc_budget_mean_zonal(contOCN,zdivOCN,hdivOCN);
    zm_heat_ice=calc_budget_mean_zonal(contICE,zdivICE,hdivICE);

    %ocean diffusion alone
    tmpU=myparms.rcp*(DFxE_TH); tmpV=myparms.rcp*(DFyE_TH);
    hdivOCN=calc_UV_conv(tmpU,tmpV);
    zm_heat_ocn_diff=calc_budget_mean_zonal(0*contOCN,0*zdivOCN,hdivOCN);
    
    %=======SALT=======
    
    contOCN=myparms.rhoconst*SALT-myparms.rhoconst*AB_gS;
    contICE=myparms.SIsal0*myparms.rhoi*SIheff;
    contTOT=contOCN+contICE;
    %vertical divergence (air-sea fluxes or vertical adv/dif)
    zdivOCN=SFLUX+oceSPflx;
    zdivICE=-zdivOCN+SRELAX;
    %in linear surface we omit :
    if ~myparms.useNLFS; zdivOCN=zdivOCN-myparms.rhoconst*WSLTMASS; end;
    %working approach for real fresh water (?) and virtual salt flux
    if ~myparms.useRFWF|~myparms.useNLFS; zdivICE=-oceSflux; end;
    %for deep ocean layer :
    if kBudget>1;
        zdivOCN=-(ADVr_SLT+DFrE_SLT+DFrI_SLT)./mygrid.RAC*myparms.rhoconst;
        zdivOCN=zdivOCN+oceSPtnd;%.*msk;
    end;
    zdivTOT=zdivOCN+zdivICE;
    %horizontal divergence (advection and diffusion)
    tmpU=myparms.rhoconst*(ADVx_SLT+DFxE_SLT); tmpV=myparms.rhoconst*(ADVy_SLT+DFyE_SLT);
    hdivOCN=calc_UV_conv(tmpU,tmpV);
    tmpU=myparms.SIsal0*(myparms.rhoi*DFxEHEFF+myparms.rhoi*ADVxHEFF);
    tmpV=myparms.SIsal0*(myparms.rhoi*DFyEHEFF+myparms.rhoi*ADVyHEFF);
    hdivICE=calc_UV_conv(tmpU,tmpV); %no dh needed here
    hdivTOT=hdivOCN+hdivICE;
    
    zm_salt_tot=calc_budget_mean_zonal(contTOT,zdivTOT,hdivTOT);
    zm_salt_ocn=calc_budget_mean_zonal(contOCN,zdivOCN,hdivOCN);
    zm_salt_ice=calc_budget_mean_zonal(contICE,zdivICE,hdivICE);

    %ocean diffusion alone
    tmpU=myparms.rhoconst*(DFxE_SLT); tmpV=myparms.rhoconst*(DFyE_SLT);
    hdivOCN=calc_UV_conv(tmpU,tmpV);
    zm_salt_ocn_diff=calc_budget_mean_zonal(0*contOCN,0*zdivOCN,hdivOCN);
    
%===================== COMPUTATIONAL SEQUENCE ENDS =========================%
%===================== PLOTTING SEQUENCE BEGINS    =========================%

elseif userStep==-1;%plotting

    if isempty(setDiagsParams);
      choicePlot={'all'};
    elseif isnumeric(setDiagsParams{1})&length(setDiagsParams)==1;
      choicePlot={'all'};
    elseif isnumeric(setDiagsParams{1});
      choicePlot={setDiagsParams{2:end}};
    else;
      choicePlot=setDiagsParams;
    end;

    tt=[1:length(alldiag.listTimes)];
    TT=alldiag.listTimes(tt);
    nt=length(TT);

    if (kBudget==1)&(sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'mass')));

        %1.1) ocean+seaice mass budgets
        %------------------------------
        figureL;
        %volume budget:
        subplot(2,1,1); disp_budget_mean_zonal(mygrid.LATS,alldiag.zm_vol_tot,'kg/m2','Mass (incl. ice)');
        %cumulative integral:
        tmp1=repmat(alldiag.zm_area,[3 1 1]);
        tmp1=tmp1.*alldiag.zm_vol_tot; cumbudg=cumsum(tmp1,2); 
        subplot(2,1,2); disp_budget_mean_zonal(mygrid.LATS,cumbudg,'kg','Mass (incl. ice)');
        %add to tex file
        myCaption={myYmeanTxt,'mass budget (ocean+ice) at each latitude in kg/m2 (upper) and integrated from South (lower).'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); end;

        %1.2) ice mass budgets
        %---------------------
        figureL;
        %volume budget:
        subplot(2,1,1); disp_budget_mean_zonal(mygrid.LATS,alldiag.zm_vol_ice,'kg/m2','Mass (only ice)');
        %cumulative integral:
        tmp1=repmat(alldiag.zm_area,[3 1 1]);
        tmp1=tmp1.*alldiag.zm_vol_ice; cumbudg=cumsum(tmp1,2);
        subplot(2,1,2); disp_budget_mean_zonal(mygrid.LATS,cumbudg,'kg','Mass (only ice)');
        %add to tex file
        myCaption={myYmeanTxt,'mass budget (only ice) at each latitude in kg/m2 (upper) and integrated from South (lower).'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); end;
    end;

    if (sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'mass')));

    %1.3) ocean mass budgets
    %-----------------------
        figureL;
        %volume budget:
        subplot(2,1,1); disp_budget_mean_zonal(mygrid.LATS,alldiag.zm_vol_ocn,'kg/m2','Mass (ocean only)');
        %cumulative integral:
        tmp1=repmat(alldiag.zm_area,[3 1 1]);
        tmp1=tmp1.*alldiag.zm_vol_ocn; cumbudg=cumsum(tmp1,2);
        subplot(2,1,2); disp_budget_mean_zonal(mygrid.LATS,cumbudg,'kg','Mass (ocean only)');
        %add to tex file
        myCaption={myYmeanTxt,'mass budget (ocean only) at each latitude in kg/m2 (upper) and integrated from South (lower).'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); end;

    end;

    if (kBudget==1)&(sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'heat')));

        %1.1) ocean+seaice heat budgets
        %------------------------------
        figureL;
        %heat budget:
        subplot(2,1,1); disp_budget_mean_zonal(mygrid.LATS,alldiag.zm_heat_tot,'J/m2','Heat (incl. ice)');
        %cumulative integral:
        tmp1=repmat(alldiag.zm_area,[3 1 1]);
        tmp1=tmp1.*alldiag.zm_heat_tot; cumbudg=cumsum(tmp1,2);
        subplot(2,1,2); disp_budget_mean_zonal(mygrid.LATS,cumbudg,'J','Heat (incl. ice)');
        %add to tex file
        myCaption={myYmeanTxt,'heat budget (ocean+ice) at each latitude in J/m2 (upper) and integrated from South (lower).'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); end;

        %1.2) ice heat budgets
        %---------------------
        figureL;
        %heat budget:
        subplot(2,1,1); disp_budget_mean_zonal(mygrid.LATS,alldiag.zm_heat_ice,'J/m2','Heat (only ice)');
        %cumulative integral:
        tmp1=repmat(alldiag.zm_area,[3 1 1]);
        tmp1=tmp1.*alldiag.zm_heat_ice; cumbudg=cumsum(tmp1,2);
        subplot(2,1,2); disp_budget_mean_zonal(mygrid.LATS,cumbudg,'J','Heat (only ice)');
        %add to tex file
        myCaption={myYmeanTxt,'heat budget (only ice) at each latitude in J/m2 (upper) and integrated from South (lower).'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); end;
    end;

    if (sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'heat')));

    %1.3) ocean heat budgets
    %-----------------------
        figureL;
        %heat budget:
        subplot(2,1,1); disp_budget_mean_zonal(mygrid.LATS,alldiag.zm_heat_ocn,'J/m2','Heat (ocean only)');
        %cumulative integral:
        tmp1=repmat(alldiag.zm_area,[3 1 1]);
        tmp1=tmp1.*alldiag.zm_heat_ocn; cumbudg=cumsum(tmp1,2);
        subplot(2,1,2); disp_budget_mean_zonal(mygrid.LATS,cumbudg,'J','Heat (ocean only)');
        %add to tex file
        myCaption={myYmeanTxt,'heat budget (ocean only) at each latitude in J/m2 (upper) and integrated from South (lower).'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); end;

    end;

    if (kBudget==1)&(sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'salt')));

        %1.1) ocean+seaice salt budgets
        %------------------------------
        figureL;
        %salt budget:
        subplot(2,1,1); disp_budget_mean_zonal(mygrid.LATS,alldiag.zm_salt_tot,'g/m2','Salt (incl. ice)');
        %cumulative integral:
        tmp1=repmat(alldiag.zm_area,[3 1 1]);
        tmp1=tmp1.*alldiag.zm_salt_tot; cumbudg=cumsum(tmp1,2);
        subplot(2,1,2); disp_budget_mean_zonal(mygrid.LATS,cumbudg,'g','Salt (incl. ice)');
        %add to tex file
        myCaption={myYmeanTxt,'salt budget (ocean+ice) at each latitude in g/m2 (upper) and integrated from South (lower).'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); end;

        %1.2) ice salt budgets
        %---------------------
        figureL;
        %salt budget:
        subplot(2,1,1); disp_budget_mean_zonal(mygrid.LATS,alldiag.zm_salt_ice,'g/m2','Salt (only ice)');
        %cumulative integral:
        tmp1=repmat(alldiag.zm_area,[3 1 1]);
        tmp1=tmp1.*alldiag.zm_salt_ice; cumbudg=cumsum(tmp1,2);
        subplot(2,1,2); disp_budget_mean_zonal(mygrid.LATS,cumbudg,'g','Salt (only ice)');
        %add to tex file
        myCaption={myYmeanTxt,'salt budget (only ice) at each latitude in g/m2 (upper) and integrated from South (lower).'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); end;
    end;

    if (sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'salt')));

    %1.3) ocean salt budgets
    %-----------------------
        figureL;
        %salt budget:
        subplot(2,1,1); disp_budget_mean_zonal(mygrid.LATS,alldiag.zm_salt_ocn,'g/m2','Salt (ocean only)');
        %cumulative integral:
        tmp1=repmat(alldiag.zm_area,[3 1 1]);
        tmp1=tmp1.*alldiag.zm_salt_ocn; cumbudg=cumsum(tmp1,2);
        subplot(2,1,2); disp_budget_mean_zonal(mygrid.LATS,cumbudg,'g','Salt (ocean only)');
        %add to tex file
        myCaption={myYmeanTxt,'salt budget (ocean only) at each latitude in g/m2 (upper) and integrated from South (lower).'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); end;

    end;

end;

