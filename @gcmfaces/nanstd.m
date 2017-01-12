function r = nanstd(p,varargin)
%overloaded gcmfaces nanstd function :
%  1) if single gcmfaces argument, then returns the global nanstd over all faces
%  2) if more than one argument, then simply calls double nanstd function for 
%     each face data, passing over the other arguments

if nargin==1;
   tmp1=[];
   for iFace=1:p.nFaces;
      iF=num2str(iFace);
      eval(['tmp1=[tmp1;p.f' iF '(:)];']);
   end;
   r=nanstd(tmp1);
   return;
end;

r=p;

for iFace=1:r.nFaces;
   iF=num2str(iFace); 
   eval(['r.f' iF '=nanstd(p.f' iF ',varargin{:});']); 
end;


