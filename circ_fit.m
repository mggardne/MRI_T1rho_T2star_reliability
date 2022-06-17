function [r,xc,yc] = circ_fit(x,y);
%CIRC_FIT Fits a circle to X and Y point coordinate data.
%
%         [R,XC,YC] = CIRC_FIT(X,Y) given the X and Y coordinates of a
%         set of points, calculates a least squares fit of a circle to
%         the points.  The radius, R, and the coordinates of the center
%         of the circle (XC, YC) are returned.
%
%         NOTES:  1.  This program assumes the points lay in a plane.
%
%                 2.  Must have at least three (3) points.
%
%                 3.  See CCircle.odt or CCircle.pdf for derivation.
%
%         27-May-2010 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error([' *** ERROR in CIRC_FIT:  The X and Y coordinates of the ', ...
         'points to be fit are required as inputs!']);
end
%
% Get Data Matrix
%
xy = [x(:) y(:)];
xym = mean(xy);
npts = size(xy,1);
if npts<3
  error(' *** ERROR in CIRC_FIT:  Not enough data points!');
end
xy = xy-repmat(xym,npts,1);            % Center data
%
% Solve for Center
%
xy1 = xy(1:npts-1,:);
xy2 = xy(2:npts,:);
lhsc = xy2-xy1;
rhsc = xy2.*xy2-xy1.*xy1;
rhsc = sum(rhsc,2)/2;
xyc = lhsc\rhsc;
%
% Get Radius
%
lhsr = ones(npts,1);
rhsr = xy-repmat(xyc',npts,1);
rhsr = rhsr.*rhsr;
rhsr = sum(rhsr,2);
r2 = lhsr\rhsr;
r = sqrt(r2);
%
% Restore Offset to Data
%
xyc = xyc+xym';                        % Add back the data offset
xc = xyc(1);
yc = xyc(2);
%
end