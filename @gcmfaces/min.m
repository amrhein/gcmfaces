function r = min(p,varargin)
%overloaded gcmfaces min function :
%  1) if single gcmfaces argument, then returns the global min over all faces
%  2) if two gcmfaces arguments, then returns the min of the two at each point
%  3) otherwise calls double min function for each face, passing over the other arguments

if nargin==1;
   tmp1=[];
   for iFace=1:p.nFaces;
      iF=num2str(iFace);
      eval(['tmp1=[tmp1;p.f' iF '(:)];']);
   end;
   r=min(tmp1);
   return;
end;

if isa(varargin{1},'gcmfaces');
   r=p;
   for iFace=1:r.nFaces;
      iF=num2str(iFace);
      eval(['r.f' iF '=min(p.f' iF ',varargin{1}.f' iF ');']);
   end;
   return;
end;

r=p;
for iFace=1:r.nFaces;
   iF=num2str(iFace); 
   eval(['r.f' iF '=min(p.f' iF ',varargin{:});']); 
end;


