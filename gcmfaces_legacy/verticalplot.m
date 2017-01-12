%22/11/2003
%auteur : Gael Forget
%version : beta
%objet : a la realsiation d'un plot xz ou yz, on utilise 
%des echelles differentes au dessus/en dessous de 500m 

%remarque : old version etait spacifique a chaque type 
%	de plot, cette nouvelle se vaut generique...

function [varargout]=verticalplot(typeplot,cur_x,cur_z,varargin);
%variables d'entree : 
%1) typeplot est le nom du type de plot (px : pcolor)
%2) cur_x et cur_z sont les chamsp x et z
%3) varargin sont les autres options a passer a la routine de plot


%nouvelles ordonnees : 
xplot2=cur_x; yplot2=myverticaldeform(cur_z);

%************************************
%gestion des vars de sortie (partie 1)
tmptxt='';
if nargout==1
tmptxt='[tmp1]=';
elseif nargout>1
tmptxt='[tmp1';
for kkk=2:nargout
tmptxt=[tmptxt ',tmp' num2str(kkk)];
end
tmptxt=[tmptxt ']='];
end
%************************************


%le plot en lui meme :
eval([tmptxt typeplot '(xplot2,yplot2,varargin{:});']);


%************************************
%gestion des vars de sortie (partie 2)
if nargout>=1
for kkk=1:nargout
eval(['varargout(kkk) = {tmp' num2str(kkk) '};'])
end
end
%************************************


%reglage des ytics :
newtick=[floor(min(min(cur_z))) ceil(max(max(cur_z)))];
if (newtick(1)<=500)&(newtick(2)>500)
tmp1=100*floor(newtick(1)/100); tmp2=ceil(newtick(2)/500)*500;
%newtick=[0:1000:6000];
%newtick=[tmp1:50:250]; newtick=[ newtick(1:end) [500:500:tmp2]];
%newtick=[tmp1:100:500]; newtick=[ newtick(1:end-1) [500:500:tmp2]];
newtick=[tmp1:100:1000]; newtick=[ newtick(1:end-1) [1000:500:6000]];
elseif (newtick(1)<=500)
tmp1=floor(newtick(1)/50)*50; tmp2=ceil(newtick(2)/50)*50;
newtick=[tmp1:50:tmp2];
else
tmp1=floor(newtick(1)/500)*500; tmp2=ceil(newtick(2)/500)*500;
newtick=[tmp1:500:tmp2];
end
%on les range dans l'ordre approprie :
newtick=flipdim(newtick,2);


%format texte :
newticklabel=[];
for kkk=1:length(newtick)
newticklabel=strvcat(newticklabel,num2str(newtick(kkk)));
end

%mise a jour des yticks :
newtick=myverticaldeform(newtick);
set(gca,'YTick',newtick);
set(gca,'YTickLabel',newticklabel);



