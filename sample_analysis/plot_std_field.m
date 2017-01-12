function []=plot_std_field(choiceV3orV4);
%object:    compute and display a standard deviation field
%inputs:    choiceV3orV4 ('v3' or 'v4') selects the sample GRID

input_list_check('plot_std_field',nargin);


%%%%%%%%%%%%%%%%%
%load parameters:
%%%%%%%%%%%%%%%%%
gcmfaces_global;
myenv.nctiles=1;
if strcmp(choiceV3orV4,'v4');
    nF=5;
    fileFormat='compact';
    dir0=myenv.gcmfaces_dir;
    dirGrid=fullfile(dir0,'../GRID/');
    myenv.nctilesdir=fullfile(dir0,'sample_input/nctiles_climatology/');
else;
    error('files are missing');
    nF=1;
    dir0=fullfile(myenv.gcmfaces_dir,'sample_input/');
    fileFormat='straight';
    dirGrid=fullfile(dir0,'grid_occa/');
    myenv.nctilesdir=fullfile(dir0,'nctiles_occa/');
end;

if ~isdir(myenv.nctilesdir);
    warning(['skipping plot_std_field (missing ' myenv.nctilesdir ')']);
    return;
end;

mygrid=[]; grid_load(dirGrid,nF,fileFormat,1);

%%%%%%%%%%%
%get field:
%%%%%%%%%%%
fld=read_nctiles([myenv.nctilesdir 'ETAN/ETAN'],'ETAN');
fld=std(fld,[],3); fld(find(fld==0))=NaN;
cc=[0:0.1:1]*0.10;

%%%%%%%%%%%%
%plot field:
%%%%%%%%%%%%
if ~myenv.lessplot;
    gcmfaces_sphere(fld,cc,[],[],1);
end;


