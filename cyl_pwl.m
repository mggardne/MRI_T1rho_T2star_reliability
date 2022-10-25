function [xyzp,xyz1i,xyz2i] = cyl_pwl(xyz1,xyz2,rad,ztrap,npts,ptol)
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
%          XYZP = CYL_PWL(XYZ1,XYZ2,RAD,ZTRAP) given a Z-value, ZTRAP,
%          traps for points too far from the origin (Z = 0) of the
%          cylinder.  Only intersection points within +/- the value of
%          ZTRAP will be returned.  The default value for ZTRAP is 20.
%
%          XYZP = CYL_PWL(XYZ1,XYZ2,RAD,ZTRAP,NPTS) given the number of
%          points to interpolate, NPTS, between the two or four
%          intersection points.  If there are four intersections, the
%          number of points is distributed equally between the two
%          pairs of intersections (i.e., half to each pair).
%
%          XYZP = CYL_PWL(XYZ1,XYZ2,RAD,ZTRAP,NPTS,PTOL) given a
%          tolerance, PTOL, checks that the distance between two
%          intersection points is within PTOL.  If not, the intersections 
%          points are discarded.  The default distance, PTOL, is 0.01.
%
%          NOTES:  1.  M-files cyl_lin3.m, cyl_pl3.m, and plane_fit.m
%                  must be in the current directory or path.
%
%                  2.  At least one line must have two intersections
%                  with the cylinder to form a polygon.
%
%          30-Sep-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check for Inputs
%
if nargin<3
  error(' *** ERROR in CYL_PWL:  Three input are required!');
end
%
if nargin<4
  ztrap = 20;
end
%
if nargin<5
  npts = 1;
end
%
if nargin<6
  ptol = 0.1;           % Less than a pixel
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
%
% First Line (Femoral Bone) Intersections
%
np1 = size(xyz1,1)-1;
n = zeros(np1,1);
xyz = cell(np1,1);
%
for l = 1:np1
   [n(l),xyz{l}] = cyl_lin3(xyz1(l,:),xyz1(l+1,:),rad);
   if n(l)>0
     if any(abs(xyz{l}(:,3))>ztrap)
       idz = abs(xyz{l}(:,3))>ztrap;
       n(l) = n(l)-sum(idz);
       xyz{l}(idz,:) = [];
     end
   end
   if sum(n)==2
     break;
   end
end
%
n1 = sum(n);
%
if n1==2
  dl = diff(cell2mat(xyz(logical(n))));
  dl = sqrt(dl*dl');
  if dl<ptol
    n = zeros(np1,1);
    n1 = 0;
  end
end
%
xyz1i = cell2mat(xyz(logical(n)));
%
% Check if Line Endpoints are Within the Cylinder
%
ie1 = [];               % Inside endpoint
end1 = true(1,2);
np1 = np1+1;
%
if n1==1
  r = norm(xyz1(1,1:2));
  if r<rad
    ie1 = 1;
    end1(1) = false;
  end
  r = norm(xyz1(np1,1:2));
  if r<rad
    ie1 = np1;
    end1(2) = false;
  end
end
%
idn = find(n);
%
if size(idn,1)<2
  if isempty(ie1)
    idn = [];
  else
    if ie1>idn
      idn = idn(1)+1:ie1;
    else
      idn = ie1:idn(1);
    end
  end
else
  idn = idn(1)+1:idn(2);
end
%
if n1==2||(n1==1&&~isempty(ie1))
%
  if n1==2
    xyzi1 = [xyz1i(1,:); xyz1(idn,:); xyz1i(2,:)];
  elseif ie1>1
    xyzi1 = [xyz1i(1,:); xyz1(idn,:)];
    n1 = 2;
  else
    xyzi1 = [xyz1(idn,:); xyz1i(1,:)];
    n1 = 2;
  end
%
end
%
% Second Line (Tibial Bone or Femoral Cartilage) Intersections
%
np2 = size(xyz2,1)-1;
n = zeros(np2,1);
xyz = cell(np2,1);
%
for l = 1:np2
   [n(l),xyz{l}] = cyl_lin3(xyz2(l,:),xyz2(l+1,:),rad);
   if n(l)>0
     if any(abs(xyz{l}(:,3))>ztrap)
       idz = abs(xyz{l}(:,3))>ztrap;
       n(l) = n(l)-sum(idz);
       xyz{l}(idz,:) = [];
     end
   end
   if sum(n)==2
     break;
   end
end
%
n2 = sum(n);
%
if n2==2
  dl = diff(cell2mat(xyz(logical(n))));
  dl = sqrt(dl*dl');
  if dl<ptol
    n = zeros(np2,1);
    n2 = 0;
  end
end
%
xyz2i = cell2mat(xyz(logical(n)));
%
% Check if Line Endpoints are Within the Cylinder
%
ie2 = [];               % Inside endpoint
end2 = true(1,2);
np2 = np2+1;
%
if n2==1
  r = norm(xyz2(1,1:2));
  if r<rad
    ie2 = 1;
    end2(1) = false;
  end
  r = norm(xyz2(np2,1:2));
  if r<rad
    ie2 = np2;
    end2(2) = false;
  end
end
%
idn = find(n);
%
if size(idn,1)<2
  if isempty(ie2)
    idn = [];
  else
    if ie2>idn
      idn = idn(1)+1:ie2;
    else
      idn = ie2:idn(1);
    end
  end
else
  idn = idn(1)+1:idn(2);
end
%
if n2==2||(n2==1&&~isempty(ie2))
%
  if n2==2
    xyzi2 = [xyz2i(1,:); xyz2(idn,:); xyz2i(2,:)];
  elseif ie2>1
    xyzi2 = [xyz2i(1,:); xyz2(idn,:)];
    n2 = 2;
  else
    xyzi2 = [xyz2(idn,:); xyz2i(1,:)];
    n2 = 2;
  end
%
% Check Directions of Lines - Should be Opposite Directions
%
  if n1==2
    df = xyzi1(end,:)-xyzi1(1,:);
    dt = xyzi2(end,:)-xyzi2(1,:);
    dchk = df*dt';      % Dot product - should be negative - opposite directions
%
    if dchk>0
      xyz2i = flipud(xyz2i);           % Flip intersection points
      xyzi2 = flipud(xyzi2);           % Flip piece-wise linear line
      end2 = fliplr(end2);
      dt = xyzi2(end,:)-xyzi2(1,:);
      dchk = df*dt';
      if dchk>=0
        error([' *** ERROR in CYL_PWL: Directions of the lines', ...
               ' are not in opposite directions!']);
      end
    end
  end
%
end
%
ends = [end1; end2];
ends = all(ends);
%
% Get Plane Normal
%
xyzs = [xyz1; xyz2];    % Combined piece-wise linear lines
[pxyz,nv] = plane_fit(xyzs(:,1),xyzs(:,2),xyzs(:,3));
%
% Interpolate Points and Combine into Polygon
%
if n1==2&&n2==2
  pts = cyl_pl3(rad,pxyz,nv,[xyzi1([1; end],:); xyzi2([1; end],:)], ...
                ends,npts);
  xyzp = [xyzi1; pts{1}; xyzi2; pts{2}];
elseif n1==2&&n2<2
  pts = cyl_pl3(rad,pxyz,nv,[xyz1i; mean(xyz2)],ends,npts);
  xyzp = [xyzi1; flipud(pts)];
elseif n1<2&&n2==2
  pts = cyl_pl3(rad,pxyz,nv,[xyz2i; mean(xyz1)],ends,npts);
  xyzp = [xyzi2; flipud(pts)];
end
%
return