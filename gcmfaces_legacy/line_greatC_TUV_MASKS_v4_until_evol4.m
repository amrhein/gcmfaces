function [LINES_out]=line_greatC_TUV_MASKS(varargin);

if nargin==1; doDisplay=varargin{1}; else; doDisplay=0; end;

global mygrid;

for iy=1:21;

switch iy;
case 1; lons=[-173 -164]; lats=[65.5 65.5]; name='Bering Strait';
case 2; lons=[-5 -5]; lats=[34 40]; name='Gibraltar';
case 3; lons=[-81 -77]; lats=[28 27]; name='Florida Strait';
case 4; lons=[-81 -79]; lats=[28 22]; name='Florida Strait W1';
case 5; lons=[-76 -76]; lats=[21 8]; name='Florida Strait S1';
case 6; lons=[-77 -77]; lats=[27 25]; name='Florida Strait E1';
case 7; lons=[-77 -77]; lats=[25 22]; name='Florida Strait E2';
case 8; lons=[-65 -50]; lats=[66 66]; name='Davis Strait';
case 9; lons=[-35 -20]; lats=[67 65]; name='Denmark Strait';
case 10; lons=[-16 -7]; lats=[65 62.5]; name='Iceland Feroe';
case 11; lons=[-6.5 -4]; lats=[62.5 57]; name='Feroe England';
case 12; lons=[-4 8]; lats=[57 62];  name='England Norway';
case 13; lons=[-68 -63]; lats=[-54 -66]; name='Drake Passage';
case 14; lons=[103 103]; lats=[-1 4]; name='Indonesia W1';
case 15; lons=[104 109]; lats=[-3 -8]; name='Indonesia W2';
case 16; lons=[113 118]; lats=[-8.5 -8.5]; name='Indonesia W3';
case 17; lons=[118 127 ]; lats=[-8.5 -15]; name='Indonesia W4';
case 18; lons=[127 127]; lats=[-25 -68]; name='Australia Antarctica';
case 19; lons=[38 46]; lats=[-10 -22]; name='Madagascar Channel';
case 20; lons=[46 46]; lats=[-22 -69]; name='Madagascar Antarctica';
case 21; lons=[20 20]; lats=[-30 -69.5]; name='South Africa Antarctica';
end;

%compute:
line_cur=line_greatC_TUV_mask(lons,lats);
line_cur.name=name;

if doDisplay;
   ii=5;
   tmp1=mygrid.hFacW{ii}(:,:,1); tmp2=line_cur.mmuIn{ii}; tmp1(~isnan(tmp2))=-1; figure; imagesc(tmp1);
   tmp1=mygrid.hFacS{ii}(:,:,1); tmp2=line_cur.mmvIn{ii}; tmp1(~isnan(tmp2))=-1; figure; imagesc(tmp1);

   [X,Y,MSK]=convert_llc2pcol(mygrid.XC,mygrid.YC,mygrid.Depth,1,1); MSK(MSK==0)=NaN;
   figure; pcolor(X,Y,MSK); axis([-180 180 -90 90]); shading flat; caxis([0 6000]);

   hold on; plot(lons(1),lats(1),'r.','MarkerSize',32);
   plot(lons(2),lats(2),'k.','MarkerSize',32);

   [X,Y,mmtIn]=convert_llc2pcol(mygrid.XC,mygrid.YC,line_cur.mmtIn,1,1);
   kk=find(mmtIn==1); plot(X(kk),Y(kk),'m.');
   %[X,Y,mmtOut]=convert_llc2pcol(mygrid.XC,mygrid.YC,line_cur.mmtOut,1,1);
   %kk=find(mmtOut==1); plot(X(kk),Y(kk),'cx');
end;

if iy==1;
   LINES_out=line_cur;
else;
   LINES_out(iy)=line_cur;
end;

end;


