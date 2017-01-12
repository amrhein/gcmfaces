%this function does a running mean along dimension dim_cur,
%averaging over indices [-nb_cur:nb_cur]
%
function [field_out]=nanrunmean_cycle(field_in,nb_cur,dim_cur);

if nb_cur~=0

size_cur=size(field_in); 

perm1to2=[1:length(size_cur)]; 
perm1to2=[dim_cur perm1to2(find(perm1to2~=dim_cur))];
perm2to1=[[1:dim_cur-1]+1 1 [dim_cur+1:length(size_cur)]]; 
size_cur2=size_cur(perm1to2);

field_in2=permute(field_in,perm1to2); 

field_mask2=1*~isnan(field_in2);

field_out2=zeros(size_cur2);
count_out2=zeros(size_cur2);
for tcur=-nb_cur:nb_cur
	tmp1=circshift(field_in2,[tcur zeros(1,length(size_cur2)-1)]);
	tmp2=find(~isnan(tmp1.*field_in2));
	field_out2(tmp2)=field_out2(tmp2)+tmp1(tmp2);
	count_out2(tmp2)=count_out2(tmp2)+1;
end

tmp1=find(count_out2>0);
field_out2(tmp1)=field_out2(tmp1)./count_out2(tmp1);
field_out2(count_out2==0)=NaN;

field_out=permute(field_out2,perm2to1);

else
field_out=field_in;
end%if nb_cur~=0


