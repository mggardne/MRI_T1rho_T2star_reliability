function [mask1,mask2,trunc] = cr_mask2f(roic,npx,dist,scal,tol,iplt)
%CR_MASK2F Creates two image masks based on the two-dimensional (2-D)
%          coordinates of two lines.  The region of interests (ROIs)
%          are defined as the space between the first cartilage line and
%          the midline between the cartilage and second bone lines and
%          the space between the midline and the bone line.  The bone
%          line is assumed to be slightly longer than the cartilage
%          line.  The bone coordinates are fit to a circle and both
%          cartilage and bone coordinates are transformed to polar
%          coordinates to check the relative lengths of the cartilage
%          and bone lines.  The cartilage line is truncated if greater
%          than 10% longer than the bone data.
%
%          For creating masks for femur joint cartilage.  The first
%          line is cartilage and the second longer line is the
%          underlying bone.  The bone coordinates are fit to a circle
%          and both cartilage and bone coordinates are converted to
%          polar coordinates.  Cartilage points beyond 10% of the bone
%          angular length (in radians) are not used in making the region
%          of interest.
%
%          [MASK1,MASK2] = CR_MASK2F(ROIC,NPX) given the two-dimensional
%          coordinates of two lines defining a region of interest in a
%          two element cell array (one element for each line), ROIC,
%          and the size of the image for the mask in array, NPX,
%          creates two logical array mask for the region of interest,
%          MASK1 and MASK2.  Note:  If NPX has only one element, the
%          image is assumed to be square (symmetrical) (NPX by NPX).
%
%          [MASK1,MASK2,TRUNC] = CR_MASK2F(ROIC,NPX) returns a logical
%          true in TRUNC if the cartilage data is truncated.
%
%          [MASK1,MASK2] = CR_MASK2F(ROIC,NPX,DIST,SCAL) midline points
%          must be within a distance DIST of the second line.  The
%          default value is Inf.  A one or two element scale SCAL is
%          used to scale the X and Y coordinates before comparison with
%          the distance DIST for midline points.  See midline.m.
%
%          NOTES:  1.  M-files cr_mask2.m, in_tri2d.m, lsect2.m,
%                  lsect2a.m, lsect3.m, lsect4.m, lsect5.m, midline.m, 
%                  mk2_tri_2d.m must be in the current directory or
%                  path.
%
%          24-Jan-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in CR_MASK2F:  Two inputs are required!');
end
%
roic = roic(:);
[nr,nc] = size(roic);
if nr~=2&&nc~=1
  error([' *** ERROR in CR_MASK2F:  Input cell array must have'
         ' two elements!']);
end
%
ndim = size(npx(:),1);
if ndim>2||ndim<1
  error([' *** ERROR in CR_MASK2F:  Incorrect number of image ', ...
         'dimensions!']);
end
%
if ndim==1
  npx = [npx; npx];
end
%
if nargin<3
  dist = Inf;            % No distance checking
end
%
if isempty(dist)
  dist = Inf;
end
%
if nargin<4
  scal = [1 1];         % No scaling
end
if isempty(scal)
  scal = [1 1];
end
scal = scal(:);
nr = size(scal,1);
if nr==1
  scal = [scal scal];   % Same scaling for X and Y
else
  scal = scal(1:2)';
end
%
if nargin<5
  tol = 0.1;
end
if isempty(tol)
  tol = 0.1;
end
%
if nargin<6
  iplt = false;
end
%
% Initialize Truncation Flag
%
trunc = false;
%
% Fit Circle to Bone Data and Convert to Polar Coordinates
%
xyb = roic{2};
[r,xc,yc] = circ_fit(xyb(:,1),xyb(:,2));
%
% Fit Circle to Bone Data and Convert to Polar Coordinates
%
nb = size(xyb,1);
xyb = xyb-repmat([xc,yc],nb,1);
[thb,rb] = cart2pol(xyb(:,1),xyb(:,2));
thb(thb<-pi/2) = thb(thb<-pi/2)+2*pi;  % Catch third quadrant values
%
xyc = roic{1};
nc = size(xyc,1);
xyc = xyc-repmat([xc,yc],nc,1);
[thc,rc] = cart2pol(xyc(:,1),xyc(:,2));
thc(thc<-pi/2) = thc(thc<-pi/2)+2*pi;  % Catch third quadrant values
%
% Get Minimum, Maximum and Range for Bone Angle (Theta)
%
thbmn = min(thb);
thbmx = max(thb);
rngb = thbmx-thbmn;
cutoff = 0.1*rngb;      % 10% of bone range
%
% Remove Cartilage Points Beyond 10% of Bone Angles
%
idmn = thc<thbmn-cutoff;
idmx = thc>thbmx+cutoff;
idc = idmn|idmx;
%
if any(idc)
  xyc = xyc(idmn);
  roic{1} = xyc;
  trunc = true;
end
%
% Call cr_mask2
%
[mask1,mask2] = cr_mask2(roic,npx,dist,scal,tol,iplt);
%
return