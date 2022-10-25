function [xyzp,xyz1i,xyz2i] = cyl_pwl(xyz1,xyz2,rad)
%CYL_PWL   Given the three-dimensional coordinates of two ordered
%          piece-wise linear lines laying in a plane and the radius of
%          a cylinder centered at the origin and aligned with the Z
%          axis, returns the coordinates of the polygon bounded by the
%          cylinder and the two lines. 
%
%          XYZP = CYL_PWL(XYZ1,XYZ2,RAD) given the three-dimensional
%          coordinates of two ordered piece-wise linear lines laying in
%          a plane, XYZ1 and XYZ2, and the radius of a cylinder
%          centered at the origin and aligned with the Z axis, RAD,
%          returns the coordinates of the polygon bounded by the
%          cylinder and the two lines, XYZP.  The X, Y, and Z
%          coordinates must be in the columns with each row
%          representing an individual point.  If the lines do not
%          intersect the cylinder, XYZP is an empty matrix.
%
%          [XYZP,XYZ1I,XYZ2I] = CYL_PWL(XYZ1,XYZ2,RAD) the intersection
%          points between the first piece-wise linear line and the
%          cylinder, XYZ1I, and the second piece-wise line and the
%          cylinder, XYZ2I, may also be returned.
%
%          NOTES:  1.  M-file cyl_lin3.m must be in the current
%                  directory or path.
%
%                  2.  Both lines must each have two intersections with
%                  the cylinder to form the polygon.  While possible to
%                  have just one intersection with one of the lines and
%                  still form a polygon (triangle), the code does not
%                  allow for this possibility.
%
%          30-Sep-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR in CYL_PWL:  Three input are required!');
end
%
% Check for 3D Row Vectors
%
if size(xyz1,2)~=3||size(xyz2,2)~=3
  error([' *** ERROR in CYL_PWL:  The two input lines must', ...
         ' have the three (3) coordinates in the columns of the', ...
         ' input arrays!']);
end
%
% Initialize Outputs
%
xyzp = [];
xyz1i = [];
xyz2i = [];
%
% First Line (Femoral Bone) Intersections
%
n1 = size(xyz1,1)-1;
n = zeros(n1,1);
xyz = cell(n1,1);
%
for l = 1:n1
   [n(l),xyz{l}] = cyl_lin3(xyz1(l,:),xyz1(l+1,:),rad);
   if sum(n)==2
     break;
   end
end
%
if sum(n)~=2
  return
end
%
xyz1i = cell2mat(xyz(logical(n)));
%
idn = find(n);
if size(idn,1)==1
  idn = [];
else
  idn = idn(1)+1:idn(2);
end
%
xyzi1 = [xyz1i(1,:); xyz1(idn,:); xyz1i(2,:)];
%
% Second Line (Tibial Bone or Femoral Cartilage) Intersections
%
n2 = size(xyz2,1)-1;
n = zeros(n2,1);
xyz = cell(n2,1);
%
for l = 1:n2
   [n(l),xyz{l}] = cyl_lin3(xyz2(l,:),xyz2(l+1,:),rad);
   if sum(n)==2
     break;
   end
end
%
if sum(n)~=2
  xyz1i = [];
  return
end
%
xyz2i = cell2mat(xyz(logical(n)));
%
idn = find(n);
if size(idn,1)==1
  idn = [];
else
  idn = idn(1)+1:idn(2);
end
%
xyzi2 = [xyz2i(1,:); xyz2(idn,:); xyz2i(2,:)];
%
% Check Directions of Lines - Should be Opposite Directions
%
df = xyzi1(end,:)-xyzi1(1,:);
dt = xyzi2(end,:)-xyzi2(1,:);
dchk = df*dt';          % Dot product - should be negative - opposite directions
%
if dchk>0
  xyzi2 = flipud(xyzi2);
  dt = xyzi2(end,:)-xyzi2(1,:);
  dchk = df*dt';
  if dchk>=0
    error([' *** ERROR in CYL_PWL: Directions of the lines are not', ...
           ' in opposite directions!']);
  end
end
%
% Create and Plot Polygon
%
xyzp = [xyzi1; xyzi2];
%
return