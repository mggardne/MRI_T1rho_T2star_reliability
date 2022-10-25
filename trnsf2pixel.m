function [pts2drb,s,r,t] = trnsf2pixel(pts3d,pts2d,pts3drb)
%TRNSF2PIXEL  Using a set of corresponding three- and two-dimensional
%             points in different coordinate systems with different
%             scales, transforms a set of three-dimensional points to
%             the two-dimensional coordinate system.
%             
%             PTS2DRB = TRNSF2PIXEL(PTS3D,PTS2D,PTS3DRB) given a set of
%             three-dimensional (3D) points, PTS3D, and corresponding
%             set of two-dimensional (2D) points, PTS2D, in different
%             coordinate systems with different scales, transforms a
%             set of 3D points, PTS3DRB, to the 2D coordinate system,
%             PTS2DRB.
%
%             [PTS2DRB,S,R,T] = TRNSF2PIXEL(PTS3D,PTS2D,PTS3DRB) returns
%             the transformation parameters scale, S, rotation matrix,
%             R, and translation vector, T. 
%             
%             NOTES:  1.  The M-file decomp.m must be in the current
%                     path or directory.
%             
%                     2.  The 2D (X and Y) and 3D (X, Y and Z)
%                     coordinates must be in columns (one point per
%                     row).
%             
%             26-Apr-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR in TRNSF2PIXEL:  Three inputs are required!');
end
%
% Get and Check the Number of 3-D and 2-D Points
%
npts = size(pts3d,1);
npts2d = size(pts2d,1);
if npts~=npts2d
  error([' *** ERROR in TRNSF2PIXEL:  Numbers of corresponding 3-D', ...
         ' and 2-D points must match!']);
end
%
% Get Transform from 3D to 2D Coordinates
%
[s,r,t] = decomp(pts3d,[pts2d zeros(npts,1)]);
%
% Transform 3D Points to 2D
%
npts3d = size(pts3drb,1);              % Number of points
pts2drb = s*pts3drb*r+repmat(t,npts3d,1);   % Transform to pixel coordinates
pts2drb = pts2drb(:,1:2);              % Convert from 3D to 2D
%
return