function [budg,area]=calc_budget_mean_zonal(cont,zconv,hconv);
%object :   compute budget of cont=zconv+hconv, on average over latitute bands
%input :    cont is the time derivative of the content (e.g. d/dt[T*drf]; per m2)
%           zconv is the convergence of vertical fluxes (e.g. TFLUX; per m2)
%           hconv is the convergence of advective and diffusive fluxes (per grid cell)
%output : budg is the vector of the 3 terms mean, over each latitude band
%         area is the vector of each latitude band
%note : for a closed budget you expect 0~(contTot-hconvTot-zconvTot)/contTot

gcmfaces_global;

warning('off','MATLAB:divideByZero');

nl=length(mygrid.LATS);
budg=zeros(3,nl);
area=zeros(1,nl);
for il=1:nl;
    if il>1&il<nl;
        tmpMin=0.5*(mygrid.LATS(il-1)+mygrid.LATS(il));
        tmpMax=0.5*(mygrid.LATS(il)+mygrid.LATS(il+1));
    elseif il==1;
        tmpMin=-Inf;
        tmpMax=0.5*(mygrid.LATS(il)+mygrid.LATS(il+1));
    elseif il==nl;
        tmpMin=0.5*(mygrid.LATS(il-1)+mygrid.LATS(il));
        tmpMax=+Inf;
    end;
    areaMask=mygrid.RAC.*mygrid.mskC(:,:,1).*(mygrid.YC>=tmpMin&mygrid.YC<tmpMax);
    %vertical integrals:
    tmpcont=nansum(cont,3)./mygrid.mskC(:,:,1);
    tmphconv=nansum(hconv,3)./(mygrid.mskC(:,:,1).*mygrid.RAC);
    tmpzconv=nansum(zconv,3)./mygrid.mskC(:,:,1);
    %horizontal integrals:
    tmpcont=nansum(tmpcont.*areaMask)/nansum(areaMask);
    tmphconv=nansum(tmphconv.*areaMask)/nansum(areaMask);
    tmpzconv=nansum(tmpzconv.*areaMask)/nansum(areaMask);
    tmparea=nansum(areaMask);
    %output result:
    if tmparea>0;
      budg(:,il)=[tmpcont;tmphconv;tmpzconv];
      area(1,il)=tmparea;
    else;
      budg(:,il)=[0;0;0];
      area(1,il)=0;
    end;
    %to check the results:
    if 0; (contTot-hconvTot-zconvTot)/contTot; end;
end;

budg(isnan(budg))=0;

warning('on','MATLAB:divideByZero');

