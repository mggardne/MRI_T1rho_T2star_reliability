function r = vec2rot3d(zvec,yvec)
%VEC2ROT3D Given a three dimensional directional vector (e.g., axis of
%          a cylinder) and an off-axis vector, computes a rotation
%          matrix to rotate points so that the directional vector is
%          aligned with the Z axis.
%
%          R = VEC2ROT3D(ZVEC,YVEC) given the three dimensional
%          directional vector, ZVEC, and the off-axis vector, YVEC,
%          computes the rotation matrix, R, which rotates points so
%          that the directional vector is aligned with the Z axis.
%
%          NOTES:  None.
%
%          27-Sep-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in VEC2ROT3D:  Two input vectors are required!');
end
%
% Check for 3D Vectors
%
zvec = zvec(:)';        % Row vector
yvec = yvec(:)';        % Row vector
%
if size(zvec,2)~=3||size(yvec,2)~=3
  error([' *** ERROR in VEC2ROT3D:  Two input vectors must be of', ...
         ' length three (3)!']);
end
%
% Normalize Vectors
%
yvec = yvec./norm(yvec);
zvec = zvec./norm(zvec);
%
% Calculate Rotation Matrix
%
xvec = cross(yvec,zvec);
xvec = xvec./norm(xvec);
yvec = cross(zvec,xvec);
yvec = yvec./norm(yvec);
%
r = [xvec; yvec; zvec]';               % Global to local rotation
%
return