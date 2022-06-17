

%
nl = size(txyz,1)-1;
dmx = 0;
%
% Loop through Line Segments
%
for k = 1:nl
   [xyzp,t] = pts2lin(txyz(k,:),txyz(k+1,:)-txyz(k,:),fxyz);
   idv = t>0&t<1;
   if any(idv)
     d = norm(xyzp(idv,:)-fxyz(idv,:));
     if d>dmx
       dmx = d;
     end
   end
end