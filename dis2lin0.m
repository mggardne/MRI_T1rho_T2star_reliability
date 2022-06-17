function [dis,idd,xyzd] = dis2lin(xyzl,xyzp,refd)
%DIS2LIN   Finds the distances from a set of points to a piece-wise
%          linear line.
%
%          DIS2LIN(XYZL,XYZP) Given a matrix with the two- or three-
%          dimensional coordinates (in columns) of the points forming a
%          piece-wise linear line, XYZL, and the coordinates (in
%          columns) of points, XYZP, finds the minimum distance from
%          the set of points to the piece-wise linear line.
%
%          DIS = DIS2LIN(XYZL,XYZP) Returns the distances in array DIS.
%
%          [DIS,IDD,XYZD] = DIS2LIN(XYZL,XYZP) Returns the logical
%          index to the points associated with the distances, IDD, and
%          the coordinates of the intersection points with the
%          piece-wise linear line, XYZD.
%
%          NOTES:  1.  The M-file pts2lin.m must be in the current path
%                  or directory.
%
%                  2.  The function works with either two- or three-
%                  dimensional coordinates.
%
%          31-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in DIS2LIN:  Two inputs are required!');
end
%
if (nargin<3)
  refd = 0;
end
%
% Initialize Arrays
%
[nl,nd] = size(xyzl);
nl = nl-1;              % Number of line segments
np = size(xyzp,1);
dis = zeros(np,1);
idd = false(np,1);
xyzd = zeros(np,nd);
%
% Loop through Line Segments
%
for k = 1:nl
   [xyzi,t] = pts2lin(xyzl(k,:),xyzl(k+1,:)-xyzl(k,:),xyzp);
   idv = t>=0&t<=1;     % Index to valid intersections within segments
   if any(idv)
     d = xyzi(idv,:)-xyzp(idv,:);
     d = sqrt(sum(d.*d,2));            % Calculate distances
     idr = d>refd;      % Find distances greater than reference distance
     d = d(idr);
     idv(idv) = idr;
     if any(idv)                       % Greater than reference distance?
       idz = dis==0;    % Saved distances equal to zeros (= 0)
       idvz = idv&idz;  % Valid distances and saved distances = 0
%
       idz = ~idz;      % Saved distances not equal to zeros (~= 0)
       idvnz = idv&idz; % Valid distances and saved distances ~= 0
       idnz = d(idvnz(idv))<dis(idvnz);
       idvnz(idvnz) = idnz;            % New distances < Saved distances ~= 0
       ida = idvz|idvnz;               % Combine zeros and nonzeros for saved distances
       idb = ida(idv);                 % Combine zeros and nonzeros for new distances
       dis(ida) = d(idb);              % New distances
%        idd(ida) = find(ida);           % Index to points with distances
       idd = idd|ida;                  % Index to points with distances
       xyzd(ida,:) = xyzi(ida,:);      % Coordinates of intersection points
       keyboard
     end
   end
end
%
return