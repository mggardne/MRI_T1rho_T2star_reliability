function [varargout] = tf_coord(rotmat,xyz0,varargin)
%TF_COORD Transforms MRI knee data cell arrays by a coordinate 
%         transformation (rotation and translation).
%
%         DATT = TF_COORD(ROTMAT,XYZ0,DAT) given the rotation matrix,
%         ROTMAT, and axes offset, XYZ0, transforms the coordinates in
%         MRI knee data cell array, DAT, by rotating and translating
%         the coordinates (txyz = xyz*rotmat+xyz0).  The transformed
%         coordinates are returned in a cell array, DATT.
%
%         [DAT1T,DAT2T,DAT3T,...] = TF_COORD(ROTMAT,XYZ0,DAT1,DAT2,DAT3,
%         ...) returns additional cell arrays that have been transformed
%         by the same coordinate transformation.
%
%         NOTES:  1.  Inputs must be a MRI knee data cell arrays with
%                 the coordinate matrices for each slice in rows in the
%                 cell array.
%
%                 2.  See coord_tf.m for a transformation that first
%                 translates and then rotates the coordinates.
%
%         18-Aug-2014 * Mack Gardner-Morse
%
%         28-Sep-2022 * Mack Gardner-Morse * Replaced loop and function
%                       call with cell array functions cellfun and
%                       mat2cell.
%

%#######################################################################
%
% Check Inputs
%
if nargin<3
  error([' *** ERROR in TF_COORD:  At least one input MRI knee', ...
         ' data cell array is required!']);
end
%
if size(rotmat,1)~=3||size(rotmat,2)~=3
  error([' *** ERROR in TF_COORD:  ROTMAT must be a 3x3 rotation', ...
         '  matrix!']);
end
%
xyz0 = xyz0(:)';
if size(xyz0,2)~=3
  error(' *** ERROR in TF_COORD:  XYZ0 must be of length three (3)!');
end
%
iloop = nargin-2;
for k = 1:iloop
   dat = varargin{k};
   if ~iscell(dat)
     error(' *** ERROR in TF_COORD:  Inputs must be cell arrays!');
   end
end
%
% Check Outputs
%
if nargout~=nargin-2
  error([' *** ERROR in TF_COORD:  Number of inputs do not match', ...
        ' the number of outputs!']);
end
%
% Transform Cell Arrays
%
for k = iloop:-1:1
%
% Get Cell Array
%
   dat = varargin{k};
%
% Get Slice Information
%
   nsd = cellfun('size',dat,1);
%
% Get and Rotate All of the Coordinates
%
   xyz = cell2mat(dat);
%
   nd = size(xyz,1);
%
   toff = repmat(xyz0,nd,1);
   txyz = xyz*rotmat+toff;
%
% Reform Cell Arrays
%
   datt = mat2cell(txyz,nsd);
%
% Output Transformed Coordinates
%   
   varargout{k} = datt;
%
end
%
return