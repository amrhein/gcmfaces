function r = median(p,varargin)
%overloaded gcmfaces median function :
%  1) if single gcmfaces argument, then returns the global median over all faces
%  2) if more than one argument, then simply calls double median function for 
%     each face data, passing over the other arguments

if nargin==1;
   tmp1=[];
   for iFace=1:p.nFaces;
      iF=num2str(iFace);
      eval(['tmp1=[tmp1;p.f' iF '(:)];']);
   end;
   r=median(tmp1);
   return;
end;

r=p;

for iFace=1:r.nFaces;
   iF=num2str(iFace); 
   eval(['r.f' iF '=median(p.f' iF ',varargin{:});']); 
end;


