function [a]=convert2vector(b);
%object:    gcmfaces to vector format conversion (if gcmfaces input)
%    or:    vector to gcmfaces format conversion (if array input)
%
%notes:     - if vector input, the gcmfaces format will be the one of mygrid.XC, so 
%           the vector input must have originally been created according to convert2vector
%           - global mygrid parameters (mygrid.XC.gridType) are used

global mygrid;

if isa(b,'gcmfaces'); do_gcmfaces2vector=1; else; do_gcmfaces2vector=0; end;

if do_gcmfaces2vector;
  bb=convert2array(b);
  a=bb(:);
else;
  bb=convert2array(mygrid.XC);
  if mod(length(b(:)),length(bb(:)))~=0;
      error('vector length is inconsistent with gcmfaces objects');
  else;
      n3=length(b(:))/length(bb(:));
  end;
  b=reshape(b,[size(bb) n3]);
  a=convert2array(b);
end;



