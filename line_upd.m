function lin = line_upd(lin,idn,idx,pts3d,pts2d)
%LINE_UPD  Updates an OsiriX segmentation CSV character line with new
%          3-D and 2-D coordinates for the segmentation.
%
%          LIN = LINE_UPD(LIN,IDN,IDX,PTS3D,PTS2D) given an OsiriX
%          segmentation CSV character line, LIN, an index to the
%          "NumOfPoints" header, IDN, an index to the commas in the
%          line, IDX, the 3-D coordinates of the segmentation, PTS3D,
%          and the 2-D pixel coordinates, PTS2D, returns the character
%          line with the new number of points and new coordinates data.
%
%          [PTS3D,PTS2D,NPTS,IDN,IDX] = PARSE_CSV(HDRS,LIN) returns
%          the number of points, NPTS, an index to the "NumOfPoints"
%          header, IDN, and an index to the commas in the line, IDX.
%
%          NOTES:  1.  The character line is expected to contain comma
%                  separated values in the format of an  OsiriX
%                  segmentation CSV file.
%
%                  2.  The X, Y and Z or X and Y coordinates must be in
%                  the columns of the point coordinates arrays.
%
%                  3.  This function does not update the image number
%                  (slice number-1), ROI measures, length, area and
%                  radius measures, ROI type, nor instance UIDs before
%                  the number of points data at the beginning of the
%                  character line.  Only the number of points and the
%                  3-D and 2-D data are updated based on the input
%                  coordinates.
%
%          29-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<5)
  error(' *** ERROR in LINE_UPD:  Five inputs are required!');
end
%
% Get and Check the Number of 3-D and 2-D Points
%
npts = size(pts3d,1);
npts2d = size(pts2d,1);
if npts~=npts2d
  error([' *** ERROR in LINE_UPD:  Numbers of 3-D and 2-D points', ...
         ' do not match!']);
end
%
% Combine 3-D and 2-D Coordinates
%
pts = [pts3d pts2d];
%
% Get First Part of the Line
%
lin1 = lin(1:idx(idn-1));              % First part of line
%
% Write Number of Points and Point Coordinates as Characters
%
lin2 = sprintf('%i',npts);
lin3 = sprintf(',%.6f,%.6f,%.6f,%.6f,%.6f',pts');
%
% Combine the Separate Parts of the Line into One Line
%
lin = [lin1 lin2 lin3];
%
return