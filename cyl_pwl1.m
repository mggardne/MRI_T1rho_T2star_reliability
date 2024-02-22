function [xyzi,xyzo,xyzc,idi,ido] = cyl_pwl1(xyz,rad,ztrap,ptol)
%CYL_PWL1  Given the three-dimensional coordinates of an ordered
%          piece-wise linear line laying in a plane and the radius of
%          a cylinder centered at the origin and aligned with the Z
%          axis, returns the segments of the piece-wise linear line
%          that are within the cylinder and the segments that are
%          outside of the cylinder.
%
%          XYZI = CYL_PWL1(XYZ,RAD) given the three-dimensional
%          coordinates of an ordered piece-wise linear line laying in a
%          plane, XYZ, and the radius of a cylinder centered at the
%          origin and aligned with the Z axis, RAD, returns the
%          coordinates of the points on the piece-wise linear line
%          within the cylinder bounded by the intersection point(s),
%          XYZI.  The X, Y, and Z coordinates must be in the columns
%          with each row representing an individual point.  If the line
%          does not intersect the cylinder, XYZI is an empty matrix.
%
%          XYZI = CYL_PWL1(XYZ,RAD,ZTRAP) given a Z-value, ZTRAP,
%          traps for points too far from the origin (Z = 0) of the
%          cylinder.  Only intersection points within +/- the value of
%          ZTRAP will be returned.  The default value for ZTRAP is 20.
%
%          XYZI = CYL_PWL1(XYZ,RAD,ZTRAP,PTOL) given a tolerance, PTOL,
%          checks that the distance between two intersection points is
%          within PTOL.  If not, the intersections points are discarded.
%          The default distance, PTOL, is 0.01.
%
%          [XYZI,XYZO] = CYL_PWL1(XYZ,RAD) returns the segment(s) of
%          the piece-wise linear line outside the cylinder including
%          the intersection point(s) in cell array, XYZO.
%
%          [XYZI,XYZO,XYZC] = CYL_PWL1(XYZ,RAD) returns the
%          intersection point(s), XYZC.
%
%          [XYZI,XYZO,XYZC,IDI,IDO] = CYL_PWL1(XYZ,RAD) returns the
%          index for the piece-wise linear line points inside the
%          cylinder, IDI, and a cell array with the index(es) for the
%          points in the segments outside of the cylinder, IDO.
%
%          NOTES:  1.  M-file cyl_lin3.m must be in the current
%                  directory or path.
%
%          20-Jul-2023 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check for Inputs
%
if nargin<2
  error(' *** ERROR in CYL_PWL1:  Two inputs are required!');
end
%
if nargin<3
  ztrap = 20;
end
%
if nargin<4
  ptol = 0.1;           % Less than a pixel
end
%
% Check for 3D Column Array
%
if size(xyz,2)~=3
  error([' *** ERROR in CYL_PWL1:  The input line must have the', ...
         ' three (3) coordinates in the columns of the input array!']);
end
%
% Initialize Outputs
%
xyzi = [];
xyzo = [];
idi = [];
ido = [];
%
% Get Line Intersections with the Cylinder by Looping Along the
% Piece-Wise Linear Line
%
np = size(xyz,1)-1;     % Number of line segments
n = zeros(np,1);
xyzic = cell(np,1);
%
for l = 1:np
   [n(l),xyzic{l}] = cyl_lin3(xyz(l,:),xyz(l+1,:),rad);
   if n(l)>0
     if any(abs(xyzic{l}(:,3))>ztrap)
       idz = abs(xyzic{l}(:,3))>ztrap;
       n(l) = n(l)-sum(idz);
       xyzic{l}(idz,:) = [];
     end
   end
   if sum(n)==2
     break;
   end
end
%
% Check Intersections
%
ni = sum(n);
%
if ni==2
  dl = diff(cell2mat(xyzic(logical(n))));
  dl = sqrt(dl*dl');
  if dl<ptol
    n = zeros(np,1);
    ni = 0;
  end
end
%
% Intersection Point(s)
%
xyzc = cell2mat(xyzic(logical(n)));    % Intersection point(s)
%
% Check if Line Endpoints are Within the Cylinder
%
if ~isempty(xyzc)
  ie1 = [];             % Inside endpoint
  np = np+1;            % Number of points
%
  if ni==1
    r = norm(xyz(1,1:2));
    if r<rad
      ie1 = 1;
    end
    r = norm(xyz(np,1:2));
    if r<rad
      ie1 = np;
    end
  end
%
% Get Piece-Wise Linear Line Points Within the Cylinder
%
  idi = find(n);
%
  if size(idi,1)<2
    if isempty(ie1)
      idi = [];
    else
      if ie1>idi
        idi = idi(1)+1:ie1;
      else
        idi = ie1:idi(1);
      end
    end
  else
    idi = idi(1)+1:idi(2);
  end
%
% Add Intersection Point(s) to Points Within the Cylinder
%
  if ni==2||(ni==1&&~isempty(ie1))
%
    if ni==2
      xyzi = [xyzc(1,:); xyz(idi,:); xyzc(2,:)];
    elseif ie1>1
      xyzi = [xyzc(1,:); xyz(idi,:)];
    else
      xyzi = [xyz(idi,:); xyzc(1,:)];
    end
%
  end
%
% Get Piece-Wise Linear Line Segments Outside of the Cylinder
%
  if any(n==2)
    idx = find(n==2);
    ido = {idx+1:np; 1:idx};
  else
    ido1 = 1:np;
    ido1 = setdiff(ido1,idi,'stable');
    idod = diff(ido1);
%
    if any(idod>1)
      idx = find(idod>1);
      ido = {ido1(idx+1:end); ido1(1:idx)};
    else
      ido = {ido1};
    end
  end
%
% Add Intersection Point(s) to Points Outside the Cylinder
%
  if size(ido,1)>1
    if ido{1}(1)==1
      xyzo = {[xyz(ido{1},:); xyzc(1,:)]; [xyzc(2,:); xyz(ido{2},:)]};
    else
      xyzo = {[xyzc(2,:); xyz(ido{1},:)]; [xyz(ido{2},:); xyzc(1,:)]};
    end
  elseif ido{1}(1)==1
    xyzo = {[xyz(ido{1},:); xyzc(1,:)]};
  else
    xyzo = {[xyzc(1,:); xyz(ido{1},:)]};
  end
%
end                     % End if intersection points
%
return