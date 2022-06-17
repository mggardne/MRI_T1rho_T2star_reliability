function [lines,hdrs,ids] = rd_csv(csvdir,csvfile,slices,irho)
%RD_CSV    Reads an OsiriX segmentation CSV file and returns the lines
%          in the file as a cell array.  Also returns the headers from
%          the first line and an index to matching slice numbers.
%
%          LINES = RD_CSV(CSVDIR,CSVFILE) given the directory name in
%          the string, CSVDIR, an OsiriX segmentation CSV file name,
%          CSVFILE, returns the lines within the file as a cell array,
%          LINES.
%
%          [LINES,HDRS] = RD_CSV(CSVDIR,CSVFILE) returns the headers,
%          HDRS, from the first line of the CSV file as a cell array.
%
%          [LINES,HDRS,IDS] = RD_CSV(CSVDIR,CSVFILE,SLICES) given slice
%          numbers, SLICES, returns the index, IDS, to the lines with
%          matching slice numbers (image number + 1).
%
%          [LINES,HDRS,IDS] = RD_CSV(CSVDIR,CSVFILE,SLICES,IRHO) given
%          the integer, IRHO, checks for "img_nums" greater than 96 and
%          subtracts one (1) from "img_nums", divides by IRHO, and adds
%          one (1).  This is to account for the digitization on the
%          first spin lock time.
%
%          NOTES:  1.  The directory path usually contains the subject
%                  number, visit number, analysis type (T1rho or T2*)
%                  and bone (femur or tibia).
%
%          28-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in RD_CSV:  Two inputs are required!');
end
%
if (nargin<3)
  slices = [];
end
%
if (nargin<4)
  irho = 1;
end
%
% Initialize Index
%
slices = slices(:);
ns = size(slices,1);
if ns>0
  ids = zeros(ns,1);
else
  ids = [];
end
%
% Open File and Read First Two Lines
%
fid = fopen(fullfile(csvdir,csvfile),'r');
lin = fgetl(fid);       % Read first line of headers
hdrs = textscan(lin,'%s','Delimiter',',');
hdrs = hdrs{1};
% idn = find(startsWith(hdrs,'NumOfPoints'));
lines = {lin};
img_nums = -1;          % Image number for header line
lin = fgetl(fid);       % Read first line of data
idl = 1;                % Index for file lines
%
% Keep Reading Data Until End of File
%
while lin~=-1
     lines = [lines; lin];
     idl = [idl; idl(end)+1];          % Index for file lines
%
     idx = strfind(lin,',');           % Get commas
%
     img_nums = [img_nums; str2double(lin(1:idx(1)-1))+1]; % Image numbers
%
%      npts = eval(lin(idx(idn-1)+1:idx(idn)-1));
%
     lin = fgetl(fid);  % Get next line
end
%
% Get Matching Image Numbers
%
if any(img_nums>96)&&irho>1
  img_nums = (img_nums-1)./irho+1;
end
%
if ns>0
  [~,idn,id] = intersect(img_nums,slices,'stable');   % Check for matching slice number
  if ~isempty(id)
    ids(id) = idl(idn);
  end
end
%
% Close File
%
fclose(fid);
%
return