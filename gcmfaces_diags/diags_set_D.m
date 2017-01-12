
%select kBudget:
if ~isempty(setDiagsParams);
  kBudget=setDiagsParams{1};
else;
  kBudget=1;
end;

doMoreBudgetOutput=0;

%override default file name:
%---------------------------
tmp1=setDiags;
if kBudget>1;
    tmp1=sprintf('D%02i',kBudget);
end;
fileMat=['diags_set_' tmp1];

if userStep==1;%diags to be computed
    listDiags=['glo_vol_ocn glo_vol_tot glo_vol_ice glo_bp'];
    listDiags=[listDiags ' north_vol_ocn north_vol_tot north_vol_ice north_bp'];
    listDiags=[listDiags ' south_vol_ocn south_vol_tot south_vol_ice south_bp'];
    listDiags=[listDiags ' glo_heat_ocn glo_heat_tot glo_heat_ice'];
    listDiags=[listDiags ' north_heat_ocn north_heat_tot north_heat_ice'];
    listDiags=[listDiags ' south_heat_ocn south_heat_tot south_heat_ice'];
    listDiags=[listDiags ' glo_salt_ocn glo_salt_tot glo_salt_ice'];
    listDiags=[listDiags ' north_salt_ocn north_salt_tot north_salt_ice'];
    listDiags=[listDiags ' south_salt_ocn south_salt_tot south_salt_ice'];

elseif userStep==2;%input files and variables
    tmp1=fullfile(dirModel,'diags',filesep,'BUDG',filesep,'budg2d_snap_set2*meta');
    test3d=isempty(dir(tmp1));
    %
    listFlds={    'ETAN','SIheff','SIhsnow','THETA   ','SALT    ','PHIBOT'};
    listFlds={listFlds{:},'SIatmFW ','oceFWflx','SItflux','TFLUX','SFLUX','oceSPflx','SRELAX'};
    listFlds={listFlds{:},'oceQnet ','SIatmQnt','SIaaflux','SIsnPrcp','SIacSubl'};
    listFlds={listFlds{:},'TRELAX','WTHMASS','WSLTMASS','oceSflux','oceQsw','oceSPtnd'};
    if kBudget>1|test3d;
        listFlds={listFlds{:},'ADVr_TH','DFrE_TH','DFrI_TH','ADVr_SLT','DFrE_SLT','DFrI_SLT','WVELMASS'};
    end;
    listFlds={listFlds{:},'SDIAG1','SDIAG2','SDIAG3'};
    listFlds={listFlds{:},'UVELMASS','VVELMASS','AB_gT','AB_gS'};
    listFlds={listFlds{:},'ADVx_TH ','ADVy_TH ','DFxE_TH ','DFyE_TH '};
    listFlds={listFlds{:},'ADVx_SLT','ADVy_SLT','DFxE_SLT','DFyE_SLT'};
    listFlds={listFlds{:},'ADVxHEFF','ADVyHEFF','DFxEHEFF','DFyEHEFF'};
    listFlds={listFlds{:},'ADVxSNOW','ADVySNOW','DFxESNOW','DFyESNOW'};
    listFldsNames=deblank(listFlds);
    %
    listFiles={'rate_budg2d_snap_set1','budg2d_hflux_set1','budg2d_zflux_set1','budg2d_zflux_set2'};
    if test3d;
        listFiles={listFiles{:},'rate_budg3d_snap_set1','budg3d_hflux_set1','budg3d_zflux_set1'};
    elseif kBudget==1;
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
    if (strcmp(setDiags,'D')&test1&test2);
        fprintf('\n abort : global and regional budgets, due to missing \n');
        fprintf(['\n   ' dirModel 'diags/BUDG/budg2d_snap_set1* \n']);
        return;
    end;

    if (strcmp(setDiags,'D')&test2);
        fprintf('\n abort : global and regional budgets, due to missing \n');
        fprintf(['\n   ' dirModel 'diags/BUDG/rate_budg2d_snap_set1* \n']);
        return;
    end;
    
    %override default file name:
    %---------------------------
    tmp1=setDiags;
    if kBudget>1;
        tmp1=sprintf('D%02i',kBudget);
    end;
    fileMat=['diags_set_' tmp1 '_' num2str(tt) '.mat'];
        
    %fill in optional fields:
    %------------------------
    if isempty(who('TRELAX')); TRELAX=0; end;
    if isempty(who('SRELAX')); SRELAX=0; end;
    if isempty(who('AB_gT')); AB_gT=0; end;
    if isempty(who('AB_gS')); AB_gS=0; end;
    if isempty(who('oceSPtnd')); oceSPtnd=0; end;
    if isempty(who('oceSPflx')); oceSPflx=0; end;
    if isempty(who('PHIBOT')); PHIBOT=0; end;

    %aliases from development phase (applies to 2012 core runs)
    %---------------------------------------------------------
    if ~isempty(who('SDIAG1')); SRELAX=SDIAG1; end;
    if ~isempty(who('SDIAG2')); SIatmFW=SDIAG2; end;
    if ~isempty(who('SDIAG3')); SItflux=SDIAG3; end;

    
    %=======MASS=========
    
    if doMoreBudgetOutput;
        %indexing and sign convention:
        %- MITgcm: fluxes are >0 downward, k=1 start at free surface
        %- here, similarly: >0 downward, k=1 free surface k=2 sea floor        
        if ~test3d;
          budgO.specs.top='free surface';
          if kBudget>1; budgO.specs.top=['interface no. ' num2str(kBudget)]; end;
          budgO.specs.bottom='sea floor';
        else;
          budgO.specs.top='interface k';
          budgO.specs.bottom='interface k+1';
        end;
        budgI.specs.top='ocn-ice to atm interface';
        budgI.specs.bottom='free surface';
        
        %here we output tendencies and fluxes in kg/s
        budgMo=budgO; budgMi=budgI;
        budgMo.specs.units='kg/s';%ocean only
        budgMi.specs.units='kg/s';%ice only
        %here we output tendencies and fluxes in Watts
        budgHo=budgO; budgHi=budgI;
        budgHo.specs.units='W';%ocean only
        budgHi.specs.units='W';%ice only
        %here we output tendencies and fluxes in g/s
        budgSo=budgO; budgSi=budgI;
        budgSo.specs.units='g/s';%ocean only
        budgSi.specs.units='g/s';%ice only        
    end;
    
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
    %
    if doMoreBudgetOutput;
        if test3d;
          tmp1=mk3D(mygrid.DRF,mygrid.hFacC).*mygrid.hFacC;
          tmp2=tmp1./mk3D(mygrid.Depth,tmp1);
          tend=tmp2.*mk3D(ETAN,tmp2)*myparms.rhoconst;
        else;
          tend=contOCN;
        end;
        budgMo.tend=mk3D(mygrid.RAC,tend).*tend;%kg/s
        budgMi.tend=mygrid.RAC.*contICE;%kg/s
    end;
        
    %vertical divergence (air-sea fluxes or vertical advection)
    zdivOCN=oceFWflx;
    zdivICE=SIatmFW-oceFWflx;
    %in virtual salt flux we omit :
    if ~myparms.useRFWF; zdivOCN=0*zdivOCN; end;
    %for deep ocean layer :
    if kBudget>1; zdivOCN=-WVELMASS*myparms.rhoconst; end;
    %
    zdivTOT=zdivOCN+zdivICE;
    %
    if doMoreBudgetOutput;
      if test3d;
        trWtop=-WVELMASS*myparms.rhoconst;
        %trWtop(:,:,1)=zdivOCN;
        trWbot=trWtop(:,:,2:length(mygrid.RC));
        trWbot(:,:,length(mygrid.RC))=0;
        %
        budgMo.trWtop=mk3D(mygrid.RAC,trWtop).*trWtop;
        budgMo.trWbot=mk3D(mygrid.RAC,trWbot).*trWbot;%kg/s
      else;
        budgMo.trWtop=mygrid.RAC.*zdivOCN; budgMo.trWbot=mygrid.RAC*0;%kg/s
      end;
      budgMi.trWtop=mygrid.RAC.*(zdivICE+zdivOCN); budgMi.trWbot=mygrid.RAC.*zdivOCN;%kg/s
    end;

    %horizontal divergence (advection and ice diffusion)
    if test3d; 
      %3D UVELMASS,VVELMASS are multiplied by DRF
      %(2D diagnostics are expectedly vertically integrated by MITgcm)
      tmp1=mk3D(mygrid.DRF,UVELMASS);
      UVELMASS=tmp1.*UVELMASS;
      VVELMASS=tmp1.*VVELMASS;
    end;
    dxg=mk3D(mygrid.DXG,VVELMASS); dyg=mk3D(mygrid.DYG,UVELMASS);
    tmpUo=myparms.rhoconst*dyg.*UVELMASS;
    tmpVo=myparms.rhoconst*dxg.*VVELMASS;
    hdivOCN=calc_UV_conv(nansum(tmpUo,3),nansum(tmpVo,3));
    tmpUi=(myparms.rhoi*DFxEHEFF+myparms.rhosn*DFxESNOW+myparms.rhoi*ADVxHEFF+myparms.rhosn*ADVxSNOW);
    tmpVi=(myparms.rhoi*DFyEHEFF+myparms.rhosn*DFyESNOW+myparms.rhoi*ADVyHEFF+myparms.rhosn*ADVySNOW);
    hdivICE=calc_UV_conv(tmpUi,tmpVi); %dh needed is alerady in DFxEHEFF etc
    hdivTOT=hdivOCN+hdivICE;
    if doMoreBudgetOutput;
        budgMo.trU=tmpUo; budgMo.trV=tmpVo;%kg/s
        budgMi.trU=tmpUi; budgMi.trV=tmpVi;%kg/s
    end;
    
    %bottom pressure for comparison:
    bp=myparms.rhoconst/9.81*PHIBOT;
    
    %compute global integrals:
    %-------------------------
    msk=mygrid.mskC(:,:,kBudget);
    glo_vol_tot=calc_budget_mean_mask(contTOT,zdivTOT,hdivTOT,msk);
    glo_vol_ocn=calc_budget_mean_mask(contOCN,zdivOCN,hdivOCN,msk);
    glo_vol_ice=calc_budget_mean_mask(contICE,zdivICE,hdivICE,msk);
    glo_bp=nansum(bp.*msk.*mygrid.RAC)/nansum(msk.*mygrid.RAC);
    
    %compute northern hemisphere integrals:
    msk=mygrid.mskC(:,:,kBudget).*(mygrid.YC>0);
    north_vol_tot=calc_budget_mean_mask(contTOT,zdivTOT,hdivTOT,msk);
    north_vol_ocn=calc_budget_mean_mask(contOCN,zdivOCN,hdivOCN,msk);
    north_vol_ice=calc_budget_mean_mask(contICE,zdivICE,hdivICE,msk);
    north_bp=nansum(bp.*msk.*mygrid.RAC)/nansum(msk.*mygrid.RAC);
    
    %and southern hemisphere integrals:
    msk=mygrid.mskC(:,:,kBudget).*(mygrid.YC<=0);
    south_vol_tot=calc_budget_mean_mask(contTOT,zdivTOT,hdivTOT,msk);
    south_vol_ocn=calc_budget_mean_mask(contOCN,zdivOCN,hdivOCN,msk);
    south_vol_ice=calc_budget_mean_mask(contICE,zdivICE,hdivICE,msk);
    south_bp=nansum(bp.*msk.*mygrid.RAC)/nansum(msk.*mygrid.RAC);
    
    %=======HEAT=======

    contOCN=myparms.rcp*THETA-myparms.rcp*AB_gT;
    contICE=-myparms.flami*(SIheff*myparms.rhoi+SIhsnow*myparms.rhosn);
    %
    if doMoreBudgetOutput;
        budgHo.tend=mk3D(mygrid.RAC,contOCN).*contOCN;%Watt
        budgHi.tend=mygrid.RAC.*contICE;%Watt
    end;
    contOCN=nansum(contOCN,3);
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
    %
    if doMoreBudgetOutput;
      if test3d;
        trWtop=-(ADVr_TH+DFrE_TH+DFrI_TH)*myparms.rcp;
        %
        dd=mygrid.RF(1:end-1);
        swfrac=0.62*exp(dd/0.6)+(1-0.62)*exp(dd/20);
        swfrac(dd<-200)=0;
        swtop=mk3D(swfrac,trWtop).*mk3D(mygrid.RAC.*oceQsw,trWtop);
        swtop(isnan(mygrid.mskC))=0;
        trWtop=trWtop+swtop;
        %
        trWtop(:,:,1)=zdivOCN.*mygrid.RAC;
        trWbot=trWtop(:,:,2:length(mygrid.RC));
        trWbot(:,:,length(mygrid.RC))=0;
        %
        budgHo.trWtop=trWtop;%Watt
        budgHo.trWbot=trWbot;%Watt
      else;
        budgHo.trWtop=mygrid.RAC.*zdivOCN; budgHo.trWbot=mygrid.RAC*0;%Watt
      end;
      budgHi.trWtop=mygrid.RAC.*(zdivICE+zdivOCN); budgHi.trWbot=mygrid.RAC.*zdivOCN;%Watt
    end;

    %horizontal divergence (advection and diffusion)
    tmpUo=myparms.rcp*(ADVx_TH+DFxE_TH); tmpVo=myparms.rcp*(ADVy_TH+DFyE_TH);
    hdivOCN=calc_UV_conv(nansum(tmpUo,3),nansum(tmpVo,3));
    tmpUi=-myparms.flami*(myparms.rhoi*DFxEHEFF+myparms.rhosn*DFxESNOW+myparms.rhoi*ADVxHEFF+myparms.rhosn*ADVxSNOW);
    tmpVi=-myparms.flami*(myparms.rhoi*DFyEHEFF+myparms.rhosn*DFyESNOW+myparms.rhoi*ADVyHEFF+myparms.rhosn*ADVySNOW);
    hdivICE=calc_UV_conv(tmpUi,tmpVi); %no dh needed here
    hdivTOT=hdivOCN+hdivICE;
    if doMoreBudgetOutput;
        budgHo.trU=tmpUo; budgHo.trV=tmpVo;%Watt
        budgHi.trU=tmpUi; budgHi.trV=tmpVi;%Watt
    end;    
    
    %compute global integrals:
    %-------------------------
    msk=mygrid.mskC(:,:,kBudget);
    glo_heat_tot=calc_budget_mean_mask(contTOT,zdivTOT,hdivTOT,msk);
    glo_heat_ocn=calc_budget_mean_mask(contOCN,zdivOCN,hdivOCN,msk);
    glo_heat_ice=calc_budget_mean_mask(contICE,zdivICE,hdivICE,msk);
    
    %compute northern hemisphere integrals:
    msk=mygrid.mskC(:,:,kBudget).*(mygrid.YC>0);
    north_heat_tot=calc_budget_mean_mask(contTOT,zdivTOT,hdivTOT,msk);
    north_heat_ocn=calc_budget_mean_mask(contOCN,zdivOCN,hdivOCN,msk);
    north_heat_ice=calc_budget_mean_mask(contICE,zdivICE,hdivICE,msk);
    
    %and southern hemisphere integrals:
    msk=mygrid.mskC(:,:,kBudget).*(mygrid.YC<=0);
    south_heat_tot=calc_budget_mean_mask(contTOT,zdivTOT,hdivTOT,msk);
    south_heat_ocn=calc_budget_mean_mask(contOCN,zdivOCN,hdivOCN,msk);
    south_heat_ice=calc_budget_mean_mask(contICE,zdivICE,hdivICE,msk);

    %=======SALT=======
        
    contOCN=myparms.rhoconst*SALT-myparms.rhoconst*AB_gS;
    contICE=myparms.SIsal0*myparms.rhoi*SIheff;
    %
    if doMoreBudgetOutput;
        budgSo.tend=mk3D(mygrid.RAC,contOCN).*contOCN;%g/s
        budgSi.tend=mygrid.RAC.*contICE;%g/s
    end;    
    contOCN=nansum(contOCN,3);
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
    %
    if doMoreBudgetOutput;
      if test3d;
        nr=length(mygrid.RC);
        trWtop=-(ADVr_SLT+DFrE_SLT+DFrI_SLT)*myparms.rhoconst;
        tmp1=mk3D(oceSPflx,oceSPtnd)-cumsum(oceSPtnd,3);
        tmp1=tmp1.*mk3D(mygrid.RAC,tmp1);
        trWtop(:,:,2:nr)=trWtop(:,:,2:nr)+tmp1(:,:,1:nr-1);
        %
        trWtop(:,:,1)=zdivOCN.*mygrid.RAC;
        trWbot=trWtop(:,:,2:length(mygrid.RC));
        trWbot(:,:,length(mygrid.RC))=0;
        %
        budgSo.trWtop=trWtop;%kg/s
        budgSo.trWbot=trWbot;%kg/s
      else;
        budgSo.trWtop=mygrid.RAC.*zdivOCN; budgSo.trWbot=mygrid.RAC*0;%kg/s
      end;
      budgSi.trWtop=0*mygrid.RAC; budgSi.trWbot=budgSo.trWtop(:,:,1);%kg/s
    end;    

    %horizontal divergence (advection and diffusion)
    tmpUo=myparms.rhoconst*(ADVx_SLT+DFxE_SLT); 
    tmpVo=myparms.rhoconst*(ADVy_SLT+DFyE_SLT);
    hdivOCN=calc_UV_conv(nansum(tmpUo,3),nansum(tmpVo,3));
    tmpUi=myparms.SIsal0*(myparms.rhoi*DFxEHEFF+myparms.rhoi*ADVxHEFF);
    tmpVi=myparms.SIsal0*(myparms.rhoi*DFyEHEFF+myparms.rhoi*ADVyHEFF);
    hdivICE=calc_UV_conv(tmpUi,tmpVi); %no dh needed here
    hdivTOT=hdivOCN+hdivICE;
    if doMoreBudgetOutput;
        budgSo.trU=tmpUo; budgSo.trV=tmpVo;%g/s
        budgSi.trU=tmpUi; budgSi.trV=tmpVi;%g/s
    end;        
    
    %compute global integrals:
    %-------------------------
    msk=mygrid.mskC(:,:,kBudget);
    glo_salt_tot=calc_budget_mean_mask(contTOT,zdivTOT,hdivTOT,msk);
    glo_salt_ocn=calc_budget_mean_mask(contOCN,zdivOCN,hdivOCN,msk);
    glo_salt_ice=calc_budget_mean_mask(contICE,zdivICE,hdivICE,msk);
    
    %compute northern hemisphere integrals:
    msk=mygrid.mskC(:,:,kBudget).*(mygrid.YC>0);
    north_salt_tot=calc_budget_mean_mask(contTOT,zdivTOT,hdivTOT,msk);
    north_salt_ocn=calc_budget_mean_mask(contOCN,zdivOCN,hdivOCN,msk);
    north_salt_ice=calc_budget_mean_mask(contICE,zdivICE,hdivICE,msk);
    
    %and southern hemisphere integrals:
    msk=mygrid.mskC(:,:,kBudget).*(mygrid.YC<=0);
    south_salt_tot=calc_budget_mean_mask(contTOT,zdivTOT,hdivTOT,msk);
    south_salt_ocn=calc_budget_mean_mask(contOCN,zdivOCN,hdivOCN,msk);
    south_salt_ice=calc_budget_mean_mask(contICE,zdivICE,hdivICE,msk);

    if doMoreBudgetOutput;
        %list of budgets to output
        listbudg={'budgMo','budgHo','budgSo'};
        if kBudget==1; listbudg={listbudg{:},'budgMi','budgHi','budgSi'}; end;
        %the actual output
        for iibudg=1:length(listbudg);
            %set directory name
            dirbudg=dirMat;
            if ~isempty(strfind(dirMat,['diags_set_' setDiags '/']))
                dirbudg=fullfile(dirMat,'..',filesep);
            end;
            sufbudg=''; 
            if kBudget>1; sufbudg=num2str(kBudget); end;
            dirbudg=fullfile(dirbudg,['diags_set_' listbudg{iibudg} sufbudg],filesep);
            %
            if ~isdir(dirbudg); mkdir(dirbudg); end;
            %set file name
            filebudg=[listbudg{iibudg} '_' num2str(tt) '.mat'];
            %output to file
            eval(['tmpbudg=' listbudg{iibudg} ';']);
            save([dirbudg filebudg],'-struct','tmpbudg');
        end;
    end;

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
        %global volume budget:
% DEA: need to replace TT with TTa?
        subplot(3,1,1); disp_budget_mean_mask(TT,alldiag.glo_vol_tot,'kg/m^2','Global Mean Mass (incl. ice)');
        %add bp:
        dt=median(diff(TT))*86400; bp=dt*cumsum(alldiag.glo_bp);
        plot(TT,bp,'k'); aa=legend; bb=get(aa,'String'); bb={bb{:},'bp'}; legend(bb,'Orientation','horizontal');
        %northern hemisphere budget:
        subplot(3,1,2); disp_budget_mean_mask(TT,alldiag.north_vol_tot,'kg/m^2','Northern Mean Mass (incl. ice)');
        %add bp:
        dt=median(diff(TT))*86400; bp=dt*cumsum(alldiag.north_bp);
        plot(TT,bp,'k'); aa=legend; bb=get(aa,'String'); bb={bb{:},'bp'}; legend(bb,'Orientation','horizontal');
        %southern hemisphere budget:
        subplot(3,1,3); disp_budget_mean_mask(TT,alldiag.south_vol_tot,'kg/m^2','Southern Mean Mass (incl. ice)');
        %add bp:
        dt=median(diff(TT))*86400; bp=dt*cumsum(alldiag.south_bp);
        plot(TT,bp,'k'); aa=legend; bb=get(aa,'String'); bb={bb{:},'bp'}; legend(bb,'Orientation','horizontal');
        %add to tex file
        myCaption={myYmeanTxt,' global (upper) north (mid) and south (lower), '};
        myCaption={myCaption{:},'mass budget (ocean+ice) in kg/m$^2$.'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); elseif ~multiTimes; close; end;

        %1.2) ice mass budgets
        %---------------------
        figureL;
        subplot(3,1,1); disp_budget_mean_mask(TT,alldiag.glo_vol_ice,'kg/m^2','Global Mean Mass (only ice)');
        dt=median(diff(TT))*86400; bp=dt*cumsum(alldiag.glo_bp);
        plot(TT,bp,'k'); aa=legend; bb=get(aa,'String'); bb={bb{:},'bp'}; legend(bb,'Orientation','horizontal');
        subplot(3,1,2); disp_budget_mean_mask(TT,alldiag.north_vol_ice,'kg/m^2','Northern Mean Mass (only ice)');
        dt=median(diff(TT))*86400; bp=dt*cumsum(alldiag.north_bp);
        plot(TT,bp,'k'); aa=legend; bb=get(aa,'String'); bb={bb{:},'bp'}; legend(bb,'Orientation','horizontal');
        subplot(3,1,3); disp_budget_mean_mask(TT,alldiag.south_vol_ice,'kg/m^2','Southern Mean Mass (only ice)');
        dt=median(diff(TT))*86400; bp=dt*cumsum(alldiag.south_bp);
        plot(TT,bp,'k'); aa=legend; bb=get(aa,'String'); bb={bb{:},'bp'}; legend(bb,'Orientation','horizontal');
        %add to tex file
        myCaption={myYmeanTxt,' global (upper) north (mid) and south (lower), '};
        myCaption={myCaption{:},'mass budget (ice only) in kg/m$^2$.'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); elseif ~multiTimes; close; end;

    end;

    if (sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'mass')));

    %1.3) ocean mass budgets
    %-----------------------
    figureL;
    %global volume budget:
    subplot(3,1,1); disp_budget_mean_mask(TT,alldiag.glo_vol_ocn,'kg/m^2','Global Mean Mass (only ocean)');
    dt=median(diff(TT))*86400; bp=dt*cumsum(alldiag.glo_bp);
    plot(TT,bp,'k'); aa=legend; bb=get(aa,'String'); bb={bb{:},'bp'}; legend(bb,'Orientation','horizontal');
    subplot(3,1,2); disp_budget_mean_mask(TT,alldiag.north_vol_ocn,'kg/m^2','Northern Mean Mass (only ocean)');
    dt=median(diff(TT))*86400; bp=dt*cumsum(alldiag.north_bp);
    plot(TT,bp,'k'); aa=legend; bb=get(aa,'String'); bb={bb{:},'bp'}; legend(bb,'Orientation','horizontal');
    subplot(3,1,3); disp_budget_mean_mask(TT,alldiag.south_vol_ocn,'kg/m^2','Southern Mean Mass (only ocean)');
    dt=median(diff(TT))*86400; bp=dt*cumsum(alldiag.south_bp);
    plot(TT,bp,'k'); aa=legend; bb=get(aa,'String'); bb={bb{:},'bp'}; legend(bb,'Orientation','horizontal');
    %add to tex file
    myCaption={myYmeanTxt,' global (upper) north (mid) and south (lower), '};
    myCaption={myCaption{:},'mass budget (ocean only) in kg/m$^2$.'};
    if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); elseif ~multiTimes; close; end;

    end;

    if (kBudget==1)&(sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'heat')));

        %2.1) ocean+seaice heat budgets
        %------------------------------
        figureL;
        subplot(3,1,1); disp_budget_mean_mask(TT,alldiag.glo_heat_tot,'J/m^2','Global Mean Ocean Heat (incl. ice)');
        subplot(3,1,2); disp_budget_mean_mask(TT,alldiag.north_heat_tot,'J/m^2','Northern Mean Ocean Heat (incl. ice)');
        subplot(3,1,3); disp_budget_mean_mask(TT,alldiag.south_heat_tot,'J/m^2','Southern Mean Ocean Heat (incl. ice)');
        %add to tex file
        myCaption={myYmeanTxt,' global (upper) north (mid) and south (lower), '};
        myCaption={myCaption{:},'heat budget (ocean+ice) in J/m$^2$.'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); elseif ~multiTimes; close; end;

        %2.2) ice heat budgets
        %---------------------
        figureL;
        subplot(3,1,1); disp_budget_mean_mask(TT,alldiag.glo_heat_ice,'J/m^2','Global Mean Ocean Heat (only ice)');
        subplot(3,1,2); disp_budget_mean_mask(TT,alldiag.north_heat_ice,'J/m^2','Northern Mean Ocean Heat (only ice)');
        subplot(3,1,3); disp_budget_mean_mask(TT,alldiag.south_heat_ice,'J/m^2','Southern Mean Ocean Heat (only ice)');
        %add to tex file
        myCaption={myYmeanTxt,' global (upper) north (mid) and south (lower), '};
        myCaption={myCaption{:},'heat budget (ice only) in J/m$^2$.'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); elseif ~multiTimes; close; end;

    end;

    if (sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'heat')));

    %2.3) ocean heat budgets
    %-----------------------
    figureL;
    subplot(3,1,1); disp_budget_mean_mask(TT,alldiag.glo_heat_ocn,'J/m^2','Global Mean Ocean Heat (only ocean)');
    subplot(3,1,2); disp_budget_mean_mask(TT,alldiag.north_heat_ocn,'J/m^2','Northern Mean Ocean Heat (only ocean)');
    subplot(3,1,3); disp_budget_mean_mask(TT,alldiag.south_heat_ocn,'J/m^2','Southern Mean Ocean Heat (only ocean)');
    %add to tex file
    myCaption={myYmeanTxt,' global (upper) north (mid) and south (lower), '};
    myCaption={myCaption{:},'heat budget (ocean only) in J/m$^2$.'};
    if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); elseif ~multiTimes; close; end;

    end;

    if (kBudget==1)&(sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'salt')));

        %3.1) ocean+seaice salt budgets
        %------------------------------
        figureL;
        subplot(3,1,1); disp_budget_mean_mask(TT,alldiag.glo_salt_tot,'g/m^2','Global Mean Ocean Salt (incl. ice)');
        subplot(3,1,2); disp_budget_mean_mask(TT,alldiag.north_salt_tot,'g/m^2','Northern Mean Ocean Salt (incl. ice)');
        subplot(3,1,3); disp_budget_mean_mask(TT,alldiag.south_salt_tot,'g/m^2','Southern Mean Ocean Salt (incl. ice)');
        %add to tex file
        myCaption={myYmeanTxt,' global (upper) north (mid) and south (lower), '};
        myCaption={myCaption{:},'salt budget (ocean+ice) in g/m$^2$.'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); elseif ~multiTimes; close; end;

        %2.2) ice salt budgets
        %---------------------
        figureL;
        subplot(3,1,1); disp_budget_mean_mask(TT,alldiag.glo_salt_ice,'g/m^2','Global Mean Ocean Salt (only ice)');
        subplot(3,1,2); disp_budget_mean_mask(TT,alldiag.north_salt_ice,'g/m^2','Northern Mean Ocean Salt (only ice)');
        subplot(3,1,3); disp_budget_mean_mask(TT,alldiag.south_salt_ice,'g/m^2','Southern Mean Ocean Salt (only ice)');
        %add to tex file
        myCaption={myYmeanTxt,' global (upper) north (mid) and south (lower), '};
        myCaption={myCaption{:},'salt budget (ice only) in g/m$^2$.'};
        if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); elseif ~multiTimes; close; end;

    end;


    if (sum(strcmp(choicePlot,'all'))|sum(strcmp(choicePlot,'salt')));

    %3.3) ocean salt budgets
    %-----------------------
    figureL;
    subplot(3,1,1); disp_budget_mean_mask(TT,alldiag.glo_salt_ocn,'g/m^2','Global Mean Ocean Salt (only ocean)');
    subplot(3,1,2); disp_budget_mean_mask(TT,alldiag.north_salt_ocn,'g/m^2','Northern Mean Ocean Salt (only ocean)');
    subplot(3,1,3); disp_budget_mean_mask(TT,alldiag.south_salt_ocn,'g/m^2','Southern Mean Ocean Salt (only ocean)');
    %add to tex file
    myCaption={myYmeanTxt,' global (upper) north (mid) and south (lower), '};
    myCaption={myCaption{:},'salt budget (ocean only) in g/m$^2$.'};
    if addToTex&multiTimes; write2tex(fileTex,2,myCaption,gcf); elseif ~multiTimes; close; end;

    end;

end;
