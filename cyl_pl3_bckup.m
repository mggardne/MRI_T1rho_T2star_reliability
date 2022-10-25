function ellpts = cyl_pl3(r,pt,nv,endpts,npts)
%CYL_PL3   Given the radius of a circular cylinder centered at the
%          origin and aligned with the Z axis, and a point in a plane
%          and the normal vector to the plane, calculates points along
%          the intersection between the cylinder and plane.
%
%          ELLPTS = CYL_PL3(R,PT,NV,ENDPTS,NPTS)  the radius of a
%          circular cylinder centered at the origin and aligned with
%          the Z axis, R, a three-dimensional (3D) point in a plane, PT,
%          the 3D normal vector to the plane, NV, three (3) or four (4)
%          points along the elliptical intersection between the plane
%          and cylinder, ENDPTS, and the number of new points to
%          generate between the end points, NPTS, returns the end
%          points with NPTS points between the end points.
%
%          If the number of ENDPTS is three (3), NPTS points are
%          returned between the first two ENDPTS on the side of the
%          elliptical intersection between the plane and cylinder in
%          the direction of the third end point.
%
%          If the number of ENDPTS is four (4), the end points must be
%          ordered and NPTS/2 points are returned between the second
%          and third end points and between the fourth and first end
%          points.
%
%          NOTES:  1.  The end points array, ENDPTS, must have a point
%                  in the row and the 3D coordinates in the columns.
%
%          03-Oct-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<5)
  error(' *** ERROR in CYL_PL3:  Five inputs are required!');
end
%
% Check for 3D Coordinates and 3D Normal Vector
%
pt = pt(:)';        % Row vector
nv = nv(:)';        % Row vector
[nendpts,ncol] = size(endpts);
%
if size(pt,2)~=3||size(nv,2)~=3||ncol~=3
  error([' *** ERROR in CYL_PL3:  The input points and normal', ...
         ' vector must be of length three (3)!']);
end
%
% Check for Three (3) or Four (4) End Points
%
if ~(nendpts==3||nendpts==4)
  error([' *** ERROR in CYL_PL3:  The number of end points must', ...
         ' be either three (3) or four (4)!']);
end
%
% Insure Number of New Points is an Integer and Get Half of the Number
%
npts = round(npts);
hnpts = round(npts/2);
%
% Normalize Normal Vector
%
nv = nv./norm(nv);
%
% Get Dot Product of the Normal Vector and Point in the Plane
%
dp = nv*pt';
%
% Get the Angle Parameter for the End Points
%
x = endpts(:,1);
y = endpts(:,2);
%
t = acos(x./r);
ts = asin(y./r);
ts = sign(ts)<0;
t(ts) = 2*pi-t(ts);     % Extend range to 2*pi (360 degrees)
%
% Calculate Points Between End Points for the Case of Three End Points
%
if nendpts==3
%
  dt0 = (t(2)-t(1))./(npts+1);
  t0 = (t(1):dt0:t(2))';
  x0 = r*cos(t0);
  y0 = r*sin(t0);
  z0 = (dp-n(1)*x0-n(2)*y0)./n(3);
%
  ellpts = [x0 y0 z0];
  pchk = ellpts(hnpts,:);
  v1 = pchk-endpts(1,:);
  v2 = endpts(3,:)-endpts(1,:);
  if v1*v2'<0           % Use the other side of the ellipse
    dt = (t(2)-sign(dt0)*2*pi-t(1))/(npts+1);
    t0 = t(1):dt:t(2);
    t0 = (t(1):dt:(t(2)-sign(dt0)*2*pi))';
    x0 = r*cos(t0);
    y0 = r*sin(t0);
    z0 = (dp-n(1)*x0-n(2)*y0)./n(3);
%
    ellpts = [x0 y0 z0];
%
  end
%
% Calculate Points Between End Points for the Case of Four End Points
%
else
%
  dt = (t(3)-t(2))./(hnpts+1);
  t1 = (t(2):dt:t(3))';
  x1 = r*cos(t1);
  y1 = r*sin(t1);
  z1 = (dp-n(1)*x1-n(2)*y1)./n(3);
%
  dt = (t(1)-t(4))./(hnpts+1);
  t2 = (t(4):dt:t(4))';
  x2 = r*cos(t2);
  y2 = r*sin(t2);
  z2 = (dp-n(1)*x2-n(2)*y2)./n(3);
%
  ellpts = [[x1 y1 z1]; [x2 y2 z2]];
%
end
%
return