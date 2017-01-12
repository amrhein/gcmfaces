
%use the field produced in sample_proccessing/example_bin_average.m
if myenv.verbose>0;
    gcmfaces_msg('* apply land mask (1/NaN) to gridded before smoothing');
end;
% fld=obsMap; fld(find(isnan(fld)))=0; fld(find(mygrid.hFacCsurf==0))=NaN;
fld=obsMap.*mygrid.mskC(:,:,1);

%%%%%%%% isotropic diffusion %%%%%%%%

%choose smoothing scale: here 3 X grid spacing
if myenv.verbose>0;
    gcmfaces_msg('* set the smoothing scale to 3 grid points');
end;
distXC=3*mygrid.DXC; distYC=3*mygrid.DYC;

%do the smoothing:
if myenv.verbose>0;
    gcmfaces_msg(['* call diffsmoth2D : apply smoothing operator ' ...
        'that consists in time stepping a diffusion equation ' ...
        'with accordingly chosen diffusivity and duration. In particular ' ...
        'diffsmoth2D illustrates gradient computations (using calc_T_grad) ' ...
        'and convergence computations (using calc_UV_conv).']);
end;
obsMap_smooth=diffsmooth2D(fld,distXC,distYC);

%display results:
if myenv.verbose>0;
    gcmfaces_msg('* crude display of results in array format');
end;

figureL; qwckplot(obsMap_smooth); caxis([-1 1]*0.4); colorbar;
title('smoothed data')

%%%%%%%% rotated diffusion %%%%%%%%

%choose smoothing scale: here 3 X grid spacing
if myenv.verbose>0;
    gcmfaces_msg('* set anisotropic and rotated smoothing');
end;
distLarge=3*sqrt(mygrid.RAC); distSmall=1*sqrt(mygrid.RAC);
fldRef=mygrid.YC;

%do the smoothing:
if myenv.verbose>0;
    gcmfaces_msg(['* call diffsmoth2Drotated : apply anisotropic smoothing operator ' ...
        'that acts preferentially along contours of a reference field (here latitude).']);
end;
obsMap_smooth2=diffsmooth2Drotated(fld,distLarge,distSmall,mygrid.YC);

%display results:
if myenv.verbose>0;
    gcmfaces_msg('* crude display of results in array format');
end;

figureL; qwckplot(obsMap_smooth2); caxis([-1 1]*0.4); colorbar;
title('zonally smoothed data')


