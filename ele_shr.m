function hp = ele_shr(econn,xyz,shr);
%ELE_SHR  Plots 3-D elements slightly shrunken for visualizing
%         individual elements.
%
%         ELE_SHR(ECONN,XYZ) Plots the elements defined by the element
%         connectivity matrix, ECONN, where each row represents one
%         element and the X, Y and Z coordinates in a three column
%         matrix, XYZ.  
%
%         HP = ELE_SHR(ECONN,XYZ) Returns a patch handle for the
%         plotted elements.
%
%         ELE_SHR(ECONN,XYZ,SHR) Optional variable, SHR, controls the
%         amount the elements are shrunk.  The variable should be
%         greater than 0 (plots just the centers of the elements) and
%         1 (elements not shrunk at all).  By default SHR = 0.75.
%         Generally, SHR should be between 0.5 and 0.85.
%
%         NOTES:  1.  Color of element patches based on the Z
%                 coordinate (height).
%
%         09-May-2016 * Mack Gardner-Morse
%

%#######################################################################
%
% Check the Inputs
%
if (nargin<2)
  error(' *** ERROR in ELE_SHR:  Not enough input arguments!');
end
%
if (nargin<3)
  shr = 0.75;
end
%
if shr<=0|shr>1
  shr = 0.75;
end
%
% Number of Elements and Vertices
%
ne = size(econn,1);     % Number of elements
nv = size(econn,2);     % Number of vertices in each element
%
% Get Shrunken Coordinates
%
xp = reshape(xyz(econn,1),ne,nv)';
yp = reshape(xyz(econn,2),ne,nv)';
zp = reshape(xyz(econn,3),ne,nv)';
xp = repmat(mean(xp),nv,1)+shr*(xp-repmat(mean(xp),nv,1));
yp = repmat(mean(yp),nv,1)+shr*(yp-repmat(mean(yp),nv,1));
zp = repmat(mean(zp),nv,1)+shr*(zp-repmat(mean(zp),nv,1));
%
% Plot elements
%
hp = patch(xp,yp,zp,zp);
view(3);
axis equal;
%
return