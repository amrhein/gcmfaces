function [budg,area]=calc_budget_mean_mask(cont,zconv,hconv,mask);
%object :   compute budget of cont=zconv+hconv, on average over masked region
%input :    cont is the time derivative of the content (e.g. d/dt[T*drf]i; per m2)
%           zconv is the convergence of vertical fluxes (e.g. TFLUX; per m2)
%           hconv is the convergence of advective and diffusive fluxes (per grid cell)
%           mask is the region of interest mask (1 to include point; 0 or NaN to exclude it)
%output : budg is the vector of the 3 terms mean, over the masked region
%         area is the vector of each latitude band
%note : for a closed budget you expect 0~(contTot-hconvTot-zconvTot)/contTot

gcmfaces_global;

mask(mask==0)=NaN;
areaMask=mygrid.RAC.*mask;

%vertical integrals:
cont=nansum(cont,3)./mask;
zconv=nansum(zconv,3)./mask;
hconv=nansum(hconv,3)./areaMask;
%horizontal integrals:
contTot=nansum(cont.*areaMask)/nansum(areaMask);
hconvTot=nansum(hconv.*areaMask)/nansum(areaMask);
zconvTot=nansum(zconv.*areaMask)/nansum(areaMask);
%output result:
budg=[contTot;hconvTot;zconvTot];
area=nansum(areaMask);
%to check the results:
if 0; (contTot-hconvTot-zconvTot)/contTot; end;

