function []=diags_driver(dirModel,dirMat,years,setDiags,doInteractive);
%object:       compute the various cost and physics
%              diagnosits from model output
%input:        dirModel is the directory containing 'diags/' or 'nctiles/'
%              dirMat is the directory where diagnostics results will be saved
%                     if isempty(dirMat) then [dirModel 'mat/'] is used by default
%              years (vector) that states years (or sets of 12 records) to compute
%                     if years=[1:4] then the first 4 years are computed
%(optional)    setDiags further specifies physical diags to be computed
%                     if setDiags='A' then one set of diags ('A') is computed
%                     By default three sets of diags ('A','B','C') are computed.
%(optional)    doInteractive=1 allows users to specify parameters interactively
%                     doInteractive = 0 (default) uses ECCO v4 parameters
%                     and omits budgets and model-data misfits analyses
%
%notes : eventually should also use dirMat, dirTex interactively for years=[];
gcmfaces_global; global myparms;

%%%%%%%%%%%%%%%
%pre-processing
%%%%%%%%%%%%%%%
myparms.yearFirst = years(1);
myparms.yearLast = years(end);

if isempty(who('doInteractive')); doInteractive=0; end;

myswitch=diags_pre_process(dirModel,dirMat,doInteractive);
% this calls diags_grid_parms.m

dirModel=[dirModel '/'];
if isempty(dirMat); dirMat=[dirModel 'mat/']; else; dirMat=[dirMat '/']; end;
%%%%%%%%%%%%%%%%%%%%
%set loop parameters
%%%%%%%%%%%%%%%%%%%%

years=years-myparms.yearFirst+1;
 if myparms.diagsAreAnnual & myparms.recInAve(2)>length(years)
    myparms.recInAve(2) = length(years);
end
if myparms.diagsAreMonthly;
    years=years(years<=myparms.recInAve(2)/12);
    lChunk=12;
elseif myparms.diagsAreAnnual;
    %    years=years(years<=myparms.recInAve(2));
    % DEA changed for decadal output
    years=years(years<=years(myparms.recInAve(2)));
    lChunk=1;
else
    error('output is not obviously in yearly or monthly format')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now do the selected computation chunk:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(who('setDiags'));
    %physical diagnostics
    % DEA changed for decadal output
    %    for myYear=years;
    for myYear=1:length(years);
        diags_select(dirModel,dirMat,setDiags,lChunk,myYear);
    end;
    
    %need profiles, cost, ctrl, budget, 'B' to work in this way too
    %need fix for : ~myswitch.doBudget and setDiags='D' or setDiags{1}='D'
    
else;
    % DEA changed for decadal output
    %    for myYear=years;
    for myYear=length(years);
        %standard physical diagnostics :
        diags_select(dirModel,dirMat,'A',lChunk,myYear);
        if myYear==years(1);
            recInAve=[myparms.recInAve(1):myparms.recInAve(2)];
            diags_select(dirModel,dirMat,'B',1,recInAve);
        end;
        diags_select(dirModel,dirMat,'C',lChunk,myYear);
        
        %budgets :
        if myswitch.doBudget;
            for kk=myparms.budgetList;
                diags_select(dirModel,dirMat,{'D',kk},lChunk,myYear);
            end;
        end;
        
        %model-data misfits :
        %  in situ profiles fit
        if myswitch.doProfiles&myYear==years(1); insitu_diags(dirMat,1); end;
        %  altimeter fit
        if myswitch.doCost&myYear==years(1); cost_altimeter(dirModel,dirMat); end;
        %  other cost terms
        if myswitch.doCost&myYear==years(1); cost_sst(dirModel,dirMat,1); end;
        if myswitch.doCost&myYear==years(1); cost_bp(dirModel,dirMat,1); end;
        if myswitch.doCost&myYear==years(1); cost_seaicearea(dirModel,dirMat,1); end;
        %  controls
        if myswitch.doCtrl&myYear==years(1); cost_xx(dirModel,dirMat,1); end;
        
    end;
    
end;

