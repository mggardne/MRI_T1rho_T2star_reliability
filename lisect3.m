function [ip,t,ierr] = lisect3(p1,v1,p2,v2,tol,clim)
%LISECT3  Finds the intersection of two closed three-dimensional (3-D)
%         lines.
%
%         IP = LISECT3(P1,V1,P2,V2) finds the intersection of two lines
%         defined by a point (P1) and a vector (V1) for the first line
%         and by a point (P2) and a vector (V2) for the second line.
%         The X, Y and Z coordinates of the intersection are returned
%         in IP.
%
%         [IP,T] = LISECT3(P1,V1,P2,V2) returns the distance along the
%         lines to the intersection.  T(1) is the normalized (0 to 1)
%         distance along the first line to the intersection.  T(2) is
%         the distance along the second line to the intersection. 
%
%         [IP,T,IERR] = LISECT3(P1,V1,P2,V2,TOL) sets IERR to true if
%         no intersection is found within tolerance (TOL).  Default
%         tolerance is 1e-8.
%
%         [IP,T,IERR] = LISECT3(P1,V1,P2,V2,TOL,CLIM) sets IERR to true
%         if the condition number of the matrix is greater than a
%         limit, CLIM.  This is usually due to parallel or nearly
%         parallel lines.  The default condition number limit is 1e+8.
%
%         NOTES:  1.  The lines are assumed to be not parallel.
%
%                 2.  All input points and vectors must be of length
%                 three (3).
%
%                 3.  See lsect.m and lsect2.m for two-dimensional
%                 line (2-D) intersections.
%
%                 4.  See lisect.m and lisect2.m for three-dimensional
%                 line (3-D) intersections with infinite lines.
%
%         07-Aug-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  error([' *** ERROR in LISECT3:  LISECT3 requires four input', ...
         ' arguments.']);
end
%
if ((nargin<5)||(isempty(tol))||(tol<=0))
  tol = 1e-8;
end
%
if ((nargin<6)||(isempty(clim))||(clim<=1))
  clim = 1e+8;
end
%
% Check Vectors
%
p1 = p1(:);
v1 = v1(:);
p2 = p2(:);
v2 = v2(:);
n1 = size(p1,1);
n2 = size(v1,1);
n3 = size(p2,1);
n4 = size(v2,1);
%
if ((n1~=3)||(n2~=3)||(n3~=3)||(n4~=3))
  error([' *** ERROR in LISECT3:  All input points and vectors' ...
         ' must be of length three (3).']);
end
%
% Set Error Code and Set Outputs to NaNs
%
ierr = true;
ip = NaN(3,1);
t = NaN(2,1);
%
% Intersection Equations
%
A = [v1(1:2) -v2(1:2)];                % Solve for intersection
b = p2(1:2)-p1(1:2);
%
if (rank(A)~=2)
  A = [v1(2:3) -v2(2:3)];
  b = p2(2:3)-p1(2:3);
end
%
if (rank(A)~=2)
  A = [v1([1 3]) -v2([1 3])];
  b = p2([1 3])-p1([1 3]);
end
%
% Check Condition Number
%
cn = cond(A);
if cn>clim
  disp([' *** WARNING in LISECT3:  Lines are parallel' ...
            ' within the condition number limit.']);
  return
end
%
% Solve for Intersection
%
tv = A\b;
%
xyz1 = tv(1)*v1+p1;                    % Coordinates of intersection
xyz2 = tv(2)*v2+p2;
%
chk = norm(xyz1-xyz2);                 % Check intersection
if (chk>tol)
  disp([' *** WARNING in LISECT3:  Lines do not intersect' ...
            ' within tolerance.']);
else
  if all(tv>=0)&&all(tv<=1)
    ip = xyz1;
    t = tv;
    ierr = false;
  end
end
%
return