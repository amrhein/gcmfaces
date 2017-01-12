function [LINES_out]=line_zonal_TUV_MASKS(Y_in);

global mygrid;

for iy=1:length(Y_in);

line_cur=line_zonal_TUV_mask(Y_in(iy));

if iy==1;
   LINES_out=line_cur;
else;
   LINES_out(iy)=line_cur;
end;

end;


