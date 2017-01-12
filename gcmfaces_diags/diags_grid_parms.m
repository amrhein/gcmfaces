function []=diags_grid_parms(listTimes,doInteractive,dirGrid);
%object :      load grid, set params, and save both to dirMat
%input :       listTimes is the time list obtained from diags_list_times
%(optional)    doInteractive=1 allows users to specify parameters interactively
%                     doInteractive = 0 (default) uses ECCO v4 parameters
%                     and omits budgets and model-data misfits analyses

%global variables
gcmfaces_global;
global myparms;

%load grid
if nargin ==2
if isempty(dir('GRID'));
    dirGrid=input('grid directory?\n');
else;
    dirGrid='GRID/';
end;
end

nF=5;
frmt='compact';
memoryLimit=0;
grid_load(dirGrid,nF,frmt,memoryLimit);

% years first and last are now defined in diags_driver using user
% input variable years
%myparms.yearFirst=2000; %first year covered by model integration
%myparms.yearLast =2060; %last year covered by model integration

%myparms.yearInAve = [(myparms.yearLast-9) myparms.yearLast];
%myparms.yearInAve=[myparms.yearFirst myparms.yearLast]; %period for time averages and variance computations

% eccov4
%myparms.timeStep =3600; %model time step for tracers
%myparms.useNLFS  =1;%2=rstar 1=nlfs 0=linear free surface
% llc45
myparms.timeStep =43200; %model time step for tracers
myparms.useNLFS  =0;%2=rstar 1=nlfs 0=linear free surface

myparms.iceModel =1;%0=use freezing point   1=use pkg/seaice   2=use pkg/thsice
myparms.useRFWF  =1;%1=real fresh water flux 0=virtual salt flux
myparms.rhoconst =1029; %sea water density
myparms.rcp      =3994*myparms.rhoconst; % sea water rho X heat capacity
myparms.rhoi     = 910; %sea ice density
myparms.rhosn    = 330; %snow density
myparms.flami    = 3.34e05; % latent heat of fusion of ice/snow (J/kg)
myparms.flamb    = 2.50e06; % latent heat of evaporation (J/kg)
myparms.SIsal0   =4;
    %myparms.yearFirst =ceil(listTimes(1)*myparms.timeStep/86400);
    %myparms.yearLast =floor(listTimes(end)*myparms.timeStep/86400); %last year covered by model integration

    % Define averaging interval. nb this must be more than a year apparently
    %if myparms.yearLast<1
    %    error('need more than one year of output')
    %elseif myparms.yearLast>9
    %    myparms.yearInAve = [(myparms.yearLast-9) myparms.yearLast];
    %else 
    %    myparms.yearInAve = [(myparms.yearFirst myparms.yearLast];
    %end

myparms.diagsNbRec=length(listTimes);
test1=median(diff(listTimes)*myparms.timeStep/86400);

% DEA changed to use on decadal avg output (treated the same as annual)
if abs(test1-30.5)<1
    myparms.diagsAreMonthly=1;
    myparms.diagsAreAnnual=0;
else
    myparms.diagsAreMonthly=0;
    myparms.diagsAreAnnual=1;
end;

%if abs(test1-30.5)<1; myparms.diagsAreMonthly=1; else; ...
%        myparms.diagsAreMonthly=0; end;
%if abs(test1-365.25)<1|abs(test1-360)<1; myparms.diagsAreAnnual=1; else; ...
%        myparms.diagsAreAnnual=0; end;

% klugy fix DA
%global FORCE_ANNUAL
%if FORCE_ANNUAL,myparms.diagsAreAnnual=0; end;

%this approximation makes things simpler:
% express listTimes (in time steps) in years:
%listTimes2=myparms.yearFirst+listTimes*myparms.timeStep/86400/365.25;
% modified to allow for pickups:
% DEA change to model year
%listTimes2=myparms.yearFirst+(listTimes-listTimes(1))*myparms.timeStep/86400/365.25;
%tmp1=-0.5*diff(listTimes,1,1)*myparms.timeStep/86400/365.25;
% This is done so that I can analyze when the first listTime is not
% 1 year (or month, or decade), e.g. from a pickup run. Note that I
% have to change the sign of tmp1 when I do this so that I am
% offsetting by a year forward instead of backwards! this ist just
% because of the difference from subtracting listTimes(1).
listTimes2=myparms.yearFirst+(listTimes-listTimes(1))*myparms.timeStep/86400/360;
tmp1=0.5*diff(listTimes,1,1)*myparms.timeStep/86400/360;
tmp1=[median(tmp1);tmp1];
listTimes2=listTimes2+tmp1;%this converts the enddate to the
                           %middate of pkg/diags

%if tmp1(end)<tmp1(end-1), tmp1(end) = []; listTimes2(end) = []; end

% define myparms.yearInAve (years in average). This should work for
% decadal and annual diags output.
%myparms.yearInAve = [listTimes2(end)-tmp1(end)-9,
%listTimes2(end)-tmp1(end)]

% ecco v4
%myparms.yearInAve = [1992,2011];

% use this one for annual!
myparms.yearInAve = [myparms.yearLast-9, myparms.yearLast]

% centennial
if median(tmp1)==50
    myparms.yearInAve = [myparms.yearLast-99,myparms.yearLast]
end
%myparms.yearInAve = [2005 2010];
%myparms.yearInAve = [2015 2020];
%myparms.yearInAve = [2001 2002];

% define myparms.recInAve (records in average)
%ii=find(listTimes2>=myparms.yearInAve(1)&listTimes2<= ...
%        myparms.yearInAve(2)+1);

% Had this before (maybe for decadal?)
%ii=find(listTimes2>=myparms.yearInAve(1)&listTimes2<= ...
%        (myparms.yearInAve(2)));

% Removed the +1 because it was extending the average period beyond
% years. This was causing problems in averaging because the cumsums
% are only computed for times in years and running a shorter
% interval after running a longer one in diags was yielding weird results!
ii=find(listTimes2>=myparms.yearInAve(1)&listTimes2<= ...
        (myparms.yearInAve(2)+1));

if myparms.diagsAreMonthly;%then restrict to full years
    ni=floor(length(ii)/12)*12; 
    if ni>0; 
      myparms.recInAve=[ii(1) ii(floor(ni))];
    else;
      myparms.recInAve=[ii(1) ii(end)];
    end;
elseif ~isempty(ii);
    myparms.recInAve=[ii(1) ii(end)];
else;
    myparms.recInAve=[1 1];
end;

