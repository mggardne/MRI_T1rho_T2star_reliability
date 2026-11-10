%
% Remove Duplicate Points and Update Connectivity
% Replace with Matlab function uniquetol?
% [xyznew, ~, ic] = uniquetol(xyz,1e-6/max(xyz(:)),'ByRows',true);
% tri_new = ic(tri);
%
  tol = 5e-4;           % 0.5 um (1/2 micron == 1/40 resolution)
  tol2 = tol*tol;
  nndup = {};
  ind = true(nnod,1);
  for k = 1:nnod
     del = xyz-repmat(xyz(k,:),nnod,1);
     del = sum(del.^2,2);
     idup = find(del<tol2);
     idup = idup(ind(idup));
     ndup = size(idup,1);
     if ndup>1          % Duplicates
       fprintf(1,'\n *** RMDUPPTS:  Duplicates found\n!');
       ind(idup(2:ndup)) = false;
       nndup = [nndup; idup'];
       for l = 2:ndup
          idx = tri==idup(l);
          tri(idx) = idup(1);
       end
     end  
  end