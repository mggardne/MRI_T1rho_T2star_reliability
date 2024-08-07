function [dmx,idd,xyzd] = mxd2linssy(xyzl,xyzp,neg)
%MXD2LINS  Finds the maximum distance from a set of points to a
%          piece-wise linear line based on the signed Y difference
%          between the lines.
%
%          MXD2LINSSY(XYZL,XYZP) Given a matrix with the two- or three-
%          dimensional coordinates (in columns) of the points forming a
%          piece-wise linear line, XYZL, and the coordinates (in
%          columns) of points, XYZP, finds the maximum distance from
%          the set of points to the piece-wise linear line.
%
%          DMX = MXD2LINSSY(XYZL,XYZP) Returns the maximum distance in
%          DMX.
%
%          DMX = MXD2LINSSY(XYZL,XYZP,NEG) Returns the maximum distance
%          with a negative Y direction if NEG is true (or nonzero)
%          (default) or the maximum distance with a positive Y
%          direction if NEG is false (or zero) in DMX.
%
%          [DMX,IDD,XYZD] = MXD2LINSSY(XYZL,XYZP) Returns the index to
%          the point with the maximum distance, IDD, and the coordinates
%          of the intersection point with the piece-wise linear line
%          from the point of maximum distance, XYZD.
%
%          NOTES:  1.  The M-file pts2lin.m must be in the current path
%                  or directory.
%
%                  2.  The function works with either two- or three-
%                  dimensional coordinates.
%
%                  3.  The function could be simplified to return all
%                  the distances from the points to the piece-wise
%                  linear line.
%
%          18-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in MXD2LINSSY:  Two inputs are required!');
end
%
if (nargin<3)||isempty(neg)
  neg = true;
end
%
% Initialize Arrays
%
nl = size(xyzl,1)-1;
dmx = 0;
idd = 0;
xyzd = zeros(1,3);
%
% Loop through Line Segments
%
for k = 1:nl
   [xyzi,t] = pts2lin(xyzl(k,:),xyzl(k+1,:)-xyzl(k,:),xyzp);
   idv = t>0&t<1;
   if any(idv)
     d = xyzi(idv,:)-xyzp(idv,:);
     ds = sign(d(:,2));
%
     if neg
       ids = ds<0;
     else
       ids = ds>0;
     end
%
     d = d(ids,:);
     idv(idv) = ids;
%
     if ~isempty(d) 
       d = sqrt(sum(d.*d,2));
       [d,idmx] = max(d);              % Find maximum
       if d>dmx                        % Update maximum?
         dmx = d;
         idd = find(idv);              % Index to point with maximum distance
         idd = idd(idmx);
         xyzd = xyzi(idd,:);           % Coordinates of intersection point
       end
     end
   end
end
%
return