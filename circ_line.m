function [n,xy] = circ_line(xy1,xy2,rad)
%CIRC_LINE Given the two end points of a two-dimensional line and the
%          radius of a circle centered at the origin, returns the
%          number of intersections and the coordinates of the
%          intersections. 
%
%          [N,XY] = CIRC_LINE(XY1,XY2,RAD) given the two end points,
%          XY1 and XY2, of a two-dimensional line and the radius, RAD,
%          of a circle centered at the origin, returns the number of
%          intersections, N (N = 0 if there are no intersections, or 1,
%          or 2 intersections), and the coordinates of the
%          intersections, XY (is empty if there are no intersections).
%          The X and Y coordinates are returned in the columns of XY.
%
%          NOTES:  1.  Also finds the intersections of a line and a
%                  cylinder centered on the origin and aligned with the
%                  Z axis.
%
%                  2.  Algorithm from:
% https://math.stackexchange.com/questions/2613781/line-cylinder-intersection?noredirect=1&lq=1
%                  See MS-Word document:  IntersectionCircleLine.docx.
%
%          29-Sep-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR in CIRC_LINE:  Three input are required!');
end
%
% Check for 2D Vectors
%
xy1 = xy1(:)';          % Row vector
xy2 = xy2(:)';          % Row vector
%
if size(xy1,2)~=2||size(xy2,2)~=2
  error([' *** ERROR in CIRC_LINE:  The two input points must', ...
         ' be of length two (2)!']);
end
%
% Calculate Intersections
%
r2 = rad*rad;           % Radius squared
%
dir = xy2-xy1;          % Line direction
dir2 = dir*dir';        % Magnitude squared
%
t = (xy1*dir')./dir2;
t = abs(t);
t2 = t*t;
%
d2 = xy1*xy1'-t2*dir2;  % Squared distance
%
% Check Intersections
%
n = 0;
xy = [];
%
if d2<=r2
%
  rk = sqrt((r2-d2)/dir2);
  tk1 = t+rk;
  tk2 = t-rk;
%
  if tk2>=0&&tk2<=1
    n = 1;
    xy = xy1+tk2*dir;
  end
%
  if tk1>=0&&tk1<=1
    if isempty(xy)
      n = 1;
      xy = xy1+tk1*dir;
    else
      n = 2;
      xy = [xy; xy1+tk1*dir];
    end
  end
%
end
%
return