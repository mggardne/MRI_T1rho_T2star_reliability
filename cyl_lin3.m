function [n,xyz] = cyl_lin3(xyz1,xyz2,rad)
%CYL_LIN3  Given the two end points of a three-dimensional line and the
%          radius of a cylinder centered at the origin and aligned with
%          the Z axis, returns the number of intersections and the
%          coordinates of the intersections. 
%
%          [N,XYZ] = CYL_LIN3(XYZ1,XYZ2,RAD) given the two end points,
%          XYZ1 and XYZ2, of a three-dimensional line and the radius,
%          RAD, of a cylinder centered at the origin and aligned with
%          the Z axis, returns the number of intersections, N (N = 0 if
%          there are no intersections, or 1, or 2 intersections), and
%          the coordinates of the intersections, XYZ (is empty if there
%          are no intersections).  The X, Y and Z coordinates are
%          returned in the columns of XYZ.
%
%          NOTES:  1.  Algorithm from:
% https://math.stackexchange.com/questions/2613781/line-cylinder-intersection?noredirect=1&lq=1
%
%          29-Sep-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR in CYL_LIN3:  Three input are required!');
end
%
% Check for 3D Row Vectors
%
xyz1 = xyz1(:)';        % Row vector
xyz2 = xyz2(:)';        % Row vector
%
if size(xyz1,2)~=3||size(xyz2,2)~=3
  error([' *** ERROR in CYL_LIN3:  The two input points must', ...
         ' be of length three (3)!']);
end
%
% Initialize Outputs
%
n = 0;
xyz = [];
%
% Get Direction Vector
%
dir = xyz2-xyz1;        % Line direction
%
% Calculate Quadratic Coefficients
%
c(3) = xyz1(1)*xyz1(1)+xyz1(2)*xyz1(2)-rad*rad;
c(2) = 2*(xyz1(1)*dir(1)+xyz1(2)*dir(2));
c(1) = dir(1)*dir(1)+dir(2)*dir(2);
%
% Get the Quadratic Roots
%
t = roots(c);
%
if ~isreal(t)
  warning('*** WARNING in CYL_LIN3:  No real roots found!');
  return
end
%
t = sort(t);            % Get line parameters in the direction of the line
%
% Check Intersections
%
if t(1)>=0&&t(1)<=1
  n = 1;
  xyz = xyz1+t(1)*dir;
end
%
if t(2)>=0&&t(2)<=1
  if isempty(xyz)
    n = 1;
    xyz = xyz1+t(2)*dir;
  else
    n = 2;
    xyz = [xyz; xyz1+t(2)*dir];
  end
end
%
return