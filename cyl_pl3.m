function ellpts = cyl_pl3(r,pt,nv,endpts,ends,npts,ang,ztrap,tol)
%CYL_PL3   Given the radius of a circular cylinder centered at the
%          origin and aligned with the Z axis, and a point in a plane
%          and the normal vector to the plane, calculates points along
%          the intersection between the cylinder and plane.
%
%          ELLPTS = CYL_PL3(R,PT,NV,ENDPTS,ENDS,NPTS) given the radius
%          of a circular cylinder centered at the origin and aligned
%          with the Z axis, R, a three-dimensional (3D) point in a
%          plane, PT, the 3D normal vector to the plane, NV, three (3)
%          or four (4) points along the elliptical intersection between
%          the plane and cylinder, ENDPTS, a two element logical array,
%          ENDS, which determines if the first and last end points
%          (element 1) and the second and third end points (element 2)
%          are to be connected by a linear line (false) or along the
%          plane and cylinder intersection (true), and the number of
%          new points to generate between the end points, NPTS, returns
%          the end points with NPTS points between the end points.
%
%          If the number of ENDPTS is three (3), NPTS points are
%          returned between the first two ENDPTS on the side of the
%          elliptical intersection between the plane and cylinder in
%          the direction of the third end point.
%
%          If the number of ENDPTS is four (4), the end points must be
%          ordered and NPTS/2 points are returned between the second
%          and third end points and between the fourth and first end
%          points.  The two sets of points are returned in a cell array.
%
%          ELLPTS = CYL_PL3(R,PT,NV,ENDPTS,NPTS,ANG) given an angle or
%          parametric distance along the parabola, ANG, checks that
%          the angle between the end points is within ANG.  If not, the
%          parametric distance is not adjusted for the full 360
%          degrees.  The default angle is 0.4189 radians (24 degrees).
%
%          ELLPTS = CYL_PL3(R,PT,NV,ENDPTS,NPTS,ANG,ZTRAP) given a
%          vertical distance along the cylinder axis (Z), ZTRAP, checks
%          that the interpolated points between the end points is
%          within ZTRAP distance of the center of the cylinder.  If not,
%          the interpolated points from the other part of the
%          intersection are used.  The default distance is 20.
%
%          ELLPTS = CYL_PL3(R,PT,NV,ENDPTS,NPTS,ANG,ZTRAP,TOL) given a
%          a tolerance, TOL, if the parametric distance between end
%          points is less than TOL, the points between the end points
%          are interpolated along a straight line instead of the
%          intersection between the plane and cylinder (conic).
%
%          NOTES:  1.  The end points array, ENDPTS, must have the
%                  points in rows and the 3D coordinates in the columns.
%
%          03-Oct-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check for Inputs
%
if nargin<6
  error(' *** ERROR in CYL_PL3:  Six inputs are required!');
end
%
if nargin<7
  ang = pi/7.5;         % 24 degrees
end
%
if nargin<8
  ztrap = 20;
end
%
if nargin<9
  tol = 5e-2;           % Much less than a pixel
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
% t90 = t>pi/2;
ts = asin(y./r);
ts2 = sign(ts)<0;
t(ts2) = 2*pi-t(ts2);   % Extend range to 2*pi (360 degrees)
%
% Calculate Points Between End Points for the Case of Three End Points
%
if nendpts==3
%
  dt = t(2)-t(1);
  if abs(dt)>tol
    dt0 = dt./(npts+1);
    t0 = (t(1):dt0:t(2))';
    x0 = r*cos(t0);
    y0 = r*sin(t0);
    z0 = (dp-nv(1)*x0-nv(2)*y0)./nv(3);
%
    ellpts = [x0 y0 z0];
%
    if any(abs(z0)>ztrap)      % Use the other side of the ellipse
      dt = (t(2)-sign(dt0)*2*pi-t(1))/(npts+1);
      t0 = (t(1):dt:(t(2)-sign(dt0)*2*pi))';
      x0 = r*cos(t0);
      y0 = r*sin(t0);
      z0 = (dp-nv(1)*x0-nv(2)*y0)./nv(3);
%
      ellpts = [x0 y0 z0];
%
    end
  else                  % Straight line
    dxyz = (endpts(2,:)-endpts(1,:))./(hnpts+1);
    x0 = (endpts(1,1):dxyz(:,1):endpts(2,1))';
    y0 = (endpts(1,2):dxyz(:,2):endpts(2,2))';
    z0 = (endpts(1,3):dxyz(:,3):endpts(2,3))';
%
    ellpts = [x0 y0 z0];
%
  end
%
  ellpts = ellpts(2:end-1,:);
%
% Calculate Points Between End Points for the Case of Four End Points
%
else
%
  dt = t(3)-t(2);
  if abs(dt)>ang
    id = [2; 3];
    t(id(ts2(id))) = t(id(ts2(id)))-2*pi;
    dt = t(3)-t(2);
  end
  if abs(dt)>tol&&ends(2)
    dt = dt./(hnpts+1);
    t1 = (t(2):dt:t(3))';
    x1 = r*cos(t1);
    y1 = r*sin(t1);
    z1 = (dp-nv(1)*x1-nv(2)*y1)./nv(3);
  else                  % Straight line
    dxyz = (endpts(3,:)-endpts(2,:))./(hnpts+1);
    x1 = (endpts(2,1):dxyz(:,1):endpts(3,1))';
    y1 = (endpts(2,2):dxyz(:,2):endpts(3,2))';
    z1 = (endpts(2,3):dxyz(:,3):endpts(3,3))';
  end
%
  dt = t(1)-t(4);
  if abs(dt)>ang
    id = [1; 4];
    t(id(ts2(id))) = t(id(ts2(id)))-2*pi;
    dt = t(1)-t(4);
  end
%
  if abs(dt)>tol&&ends(1)
    dt = dt./(hnpts+1);
    t2 = (t(4):dt:t(1))';
    x2 = r*cos(t2);
    y2 = r*sin(t2);
    z2 = (dp-nv(1)*x2-nv(2)*y2)./nv(3);
  else                  % Straight line
    dxyz = (endpts(1,:)-endpts(4,:))./(hnpts+1);
    x2 = (endpts(4,1):dxyz(:,1):endpts(1,1))';
    y2 = (endpts(4,2):dxyz(:,2):endpts(1,2))';
    z2 = (endpts(4,3):dxyz(:,3):endpts(1,3))';
  end
%
  xyz1 = [x1 y1 z1];
  xyz1 = xyz1(2:end-1,:);
%
  xyz2 = [x2 y2 z2];
  xyz2 = xyz2(2:end-1,:);
%
  ellpts = {xyz1; xyz2};
%
end
%
return