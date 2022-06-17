function [pts3d,pts2d,npts,idn,idx] = parse_csv(hdrs,lin)
%PARSE_CSV Reads an OsiriX segmentation CSV character line and returns
%          the 3-D and 2-D coordinates of the segmentation.
%
%          PTS3D = PARSE_CSV(HDRS,LIN) given a cell array of the
%          headers from the first line of the CSV file, HDRS, and a
%          character line from the CSV file, LIN, returns the 3-D
%          coordinates of the segmentation, PTS3D.  The X, Y, and Z
%          coordinates are in the three columns of the array.
%
%          [PTS3D,PTS2D] = PARSE_CSV(HDRS,LIN) returns the 2-D pixel
%          coordinates, PTS2D.  The X and Y coordinates are in the two
%          columns of the array.
%
%          [PTS3D,PTS2D,NPTS,IDN,IDX] = PARSE_CSV(HDRS,LIN) returns
%          the number of points, NPTS, an index to the "NumOfPoints"
%          header, IDN, and an index to the commas in the line, IDX.
%
%          NOTES:  1.  The character line is expected to contain comma
%                  separated values in the format of an  OsiriX
%                  segmentation CSV file.
%
%          29-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in PARSE_CSV:  Two inputs are required!');
end
%
% Get Index to "NumOfPoints" Header
%
idn = find(startsWith(hdrs,'NumOfPoints'));
%
% Get Commas in Line
%
idx = strfind(lin,',');
idx = [idx length(lin)+1];
%
% Number of Data Points for this ROI
%
npts = eval(lin(idx(idn-1)+1:idx(idn)-1));
idp = idn+(0:5:(npts-1)*5)';
pts2d = zeros(npts,2);
pts3d = zeros(npts,3);
%
% Get Point Data
%
for k = 1:npts
   pts2d(k,:) = eval(['[' lin(idx(idp(k)+3)+1:idx(idp(k)+5)-1) ']']);
   pts3d(k,:) = eval(['[' lin(idx(idp(k))+1:idx(idp(k)+3)-1) ']']);
end
%
return