function [varout]=myverticaldeform(varin)
%varin : profondeur en m (positif)
%varout : z correspondant (negatif), avec 
%	deformation selon le niveau

%var de sortie :
varout=zeros(size(varin));
%valeurs limite :
%val_lim=500;fprintf('deform v:500m\n');
val_lim=500;
%val_lim=250; fprintf('deform v:250m\n');
%val_lim=-0.0001;
val_lim2=2000;

%valeurs inferieures a 500m :
tmp1=find(varin<=val_lim);
tmp2=-100*varin(tmp1)/val_lim;
varout(tmp1)=tmp2;

%valeurs superieurs a 500m :
tmp1=find(varin>val_lim);
tmp2=-100-200*(varin(tmp1)-val_lim)./(val_lim2-val_lim);
varout(tmp1)=tmp2;



