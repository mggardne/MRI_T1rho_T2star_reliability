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
%          [DIS,IDD,XYZD] = DIS2LIN(XYZL,XYZP,REFD) Only finds distances
%          for points greater than a reference distance, REFD, from the
%          lines.  By default, the reference distance is eps^(2/3).
%
%          NOTES:  1.  The M-file pts2lin.m must be in the current path
%                  or directory.
%
%                  2.  The function works with either two- or three-
%                  dimensional coordinates.
%
%                  3.  Distances from points on the lines or within
%                  REFD distance of the lines are not returned.
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
  refd = eps^(2/3);
end
%
% Initialize Arrays
%
[nl,nd] = size(xyzl);
nl = nl-1;              % Number of line segments
np = size(xyzp,1);      % Number of points
dis = Inf(np,1);        % Valid shortest distances
idd = false(np,1);      % Index to points with valid shortest distances
xyzd = zeros(np,nd);    % Intersection points on lines to points
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
       idn = d<dis(idv);               % Find shortest distances
       idv(idv) = idn;
       dis(idv) = d(idn);              % New distances
       idd = idd|idv;                  % Index to points with distances
       xyzd(idv,:) = xyzi(idv,:);      % Coordinates of intersection points
     end
   end
end
%
dis = dis(idd);
xyzd = xyzd(idd,:);
%
return