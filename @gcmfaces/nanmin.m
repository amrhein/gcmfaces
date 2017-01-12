function r = nanmin(p,varargin)
%overloaded gcmfaces nanmin function :
%  1) if single gcmfaces argument, then returns the global nanmin over all faces
%  2) if two gcmfaces arguments, then returns the nanmin of the two at each point
%  3) otherwise calls double nanmin function for each face, passing over the other arguments

if nargin==1;
   tmp1=[];
   for iFace=1:p.nFaces;
      iF=num2str(iFace);
      eval(['tmp1=[tmp1;p.f' iF '(:)];']);
   end;
   r=nanmin(tmp1);
   return;
end;

if isa(varargin{1},'gcmfaces');
   r=p;
   for iFace=1:r.nFaces;
      iF=num2str(iFace);
      eval(['r.f' iF '=nanmin(p.f' iF ',varargin{1}.f' iF ');']);
   end;
   return;
end;

r=p;
for iFace=1:r.nFaces;
   iF=num2str(iFace); 
   eval(['r.f' iF '=nanmin(p.f' iF ',varargin{:});']); 
end;


