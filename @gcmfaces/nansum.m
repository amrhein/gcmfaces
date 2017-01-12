function r = nannansum(p,varargin)
%overloaded gcmfaces nansum function :
%  1) if single gcmfaces argument, then returns the global nansum over all faces
%  2) if more than one argument, then simply calls double nansum function for 
%     each face data, passing over the other arguments

if nargin==1;
   tmp1=[]; 
   for iFace=1:p.nFaces; 
      iF=num2str(iFace);
      eval(['tmp1=[tmp1;p.f' iF '(:)];']);
   end;
   r=nansum(tmp1);
   return;
end;

r=p;

for iFace=1:r.nFaces;
   iF=num2str(iFace); 
   eval(['r.f' iF '=nansum(p.f' iF ',varargin{:});']); 
end;


