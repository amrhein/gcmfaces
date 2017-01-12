
%test:
%lons=[-81 -77]; lats=[28 27]; name='Florida Strait';
%lons=[-81 -79]; lats=[28 22]; name='Florida Strait W1';
%lons=[-76 -76]; lats=[21 8]; name='Florida Strait S1';
%lons=[-77 -77]; lats=[27 25]; name='Florida Strait E1';
%lons=[-77 -77]; lats=[25 22]; name='Florida Strait E2';

%lons=[-65 -50]; lats=[66 66]; name='Davis Strait';
%lons=[-35 -20]; lats=[67 65]; name='Denmark Strait';
%lons=[-16 -7]; lats=[65 62.5]; name='Iceland Feroe';
%lons=[-6.5 -4]; lats=[62.5 57]; name='Feroe England';
%lons=[-4 8]; lats=[57 62];  name='England Norway';

%lons=[-68 -63]; lats=[-54 -66]; name='Drake Passage';
%lons=[103 103]; lats=[-1 4]; name='Indonesia W1';
%lons=[104 109]; lats=[-3 -8]; name='Indonesia W2';
%lons=[113 118]; lats=[-8.5 -8.5]; name='Indonesia W3';
%lons=[118 127 ]; lats=[-8.5 -15]; name='Indonesia W4';
%lons=[127 127]; lats=[-25 -68]; name='Australia Antarctica';
%lons=[38 46]; lats=[-10 -22]; name='Madagascar Channel';
%lons=[46 46]; lats=[-22 -69]; name='Madagascar Antarctica';
%lons=[20 20]; lats=[-30 -69.5]; name='South Africa Antarctica';

%compute:
line_cur=line_greatC_TUV_mask(lons,lats);
line_cur.name=name;

%check:
ii=1; close all;
tmp1=mygrid.hFacW{ii}(:,:,1); tmp2=line_cur.mmuIn{ii}; tmp1(~isnan(tmp2))=-1; figure; imagesc(tmp1);
%axis([40 60 50 100]);
%axis([50 190 40 90]);
axis([230 270 0 50]);
drawnow; refresh; pause;

tmp1=mygrid.hFacS{ii}(:,:,1); tmp2=line_cur.mmvIn{ii}; tmp1(~isnan(tmp2))=-1; figure; imagesc(tmp1);
%axis([40 60 50 100]);
%axis([50 190 40 90]);
axis([230 270 0 50]);
drawnow; refresh; pause;

[X,Y,MSK]=convert_llc2pcol(mygrid.XC,mygrid.YC,mygrid.Depth,1,1); MSK(MSK==0)=NaN;
figure; pcolor(X,Y,MSK); axis([-180 180 -90 90]); shading flat; caxis([0 6000]);
hold on; plot(lons(1),lats(1),'r.','MarkerSize',32);
plot(lons(2),lats(2),'k.','MarkerSize',32);
[X,Y,mmtIn]=convert_llc2pcol(mygrid.XC,mygrid.YC,line_cur.mmtIn,1,1);
kk=find(mmtIn==1); plot(X(kk),Y(kk),'m.');
%axis([-90 -70 10 30]);
%axis([90 140 -70 30]);
%axis([15 30 -72 -30]);
axis([-80 20 50 70]);
drawnow; refresh; pause;

