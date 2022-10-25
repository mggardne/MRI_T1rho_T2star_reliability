function [pts3drb,s,r,t] = trnsf2mm(pts2d,pts3d,pts2drb)
%TRNSF2MM  Using a set of corresponding three- and two-dimensional
%          points in different coordinate systems with different
%          scales, transforms a set of two-dimensional points to
%          the three-dimensional coordinate system.
%
%          PTS3DRB = TRNSF2MM(PTS2D,PTS3D,PTS2DRB) given a set of
%          two-dimensional (2D) points, PTS2D, and corresponding
%          set of three-dimensional (3D) points, PTS3D, in different
%          coordinate systems with different scales, transforms a
%          set of 2D points, PTS2DRB, to the 3D coordinate system,
%          PTS3DRB.
%
%          [PTS3DRB,S,R,T] = TRNSF2MM(PTS2D,PTS3D,PTS2DRB) returns the
%          transformation parameters, scale, S, rotation matrix, R, and
%          translation vector, T. 
%          
%          NOTES:  1.  The M-file decomp.m must be in the current
%                  path or directory.
%          
%                  2.  The 2D (X and Y) and 3D (X, Y and Z)
%                  coordinates must be in columns (one point per
%                  row).
%          
%          28-Apr-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR in TRNSF2MM:  Three inputs are required!');
end
%
% Get and Check the Number of 3-D and 2-D Points
%
npts2d = size(pts2d,1);
npts3d = size(pts3d,1);
if npts2d~=npts3d
  error([' *** ERROR in TRNSF2MM:  Numbers of corresponding 2-D', ...
         ' and 3-D points must match!']);
end
%
% Get Transform from 2D to 3D Coordinates
%
[s,r,t] = decomp([pts2d zeros(npts2d,1)],pts3d);
%
% Transform 2D Points to 3D
%
npts = size(pts2drb,1);
pts3drb = s*[pts2drb zeros(npts,1)]*r+repmat(t,npts,1);
%
return