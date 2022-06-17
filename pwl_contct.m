function [nid,xyzc] = pwl_contct(xyzf,xyzt,thresh)
%PWL_CONTCT Finds points within a threshold of a piece-wise linear line.
%          Also calculates the coordinates of the points on a piece-wise
%          linear line that are within a threshold of the set of points.
%
%          NID = PWL_CONTCT(XYZF,XYZT,THRESH) Given a two- or three-
%          dimensional array of coordinates of points forming a
%          piece-wise linear line, XYZF, a two- or three-dimensional
%          array of coordinates of points, XYZT, and a threshold for
%          the distance of the points from the piece-wise linear line
%          to be within, THRESHOLD, returns the index, NID, to the
%          points that are within the threshold distance to the
%          piece-wise linear line.
%
%          [NID,XYZC] = PWL_CONTCT(XYZF,XYZT,THRESH) Returns the
%          coordinates of the points on the piece-wise linear line that
%          are within a threshold of the set of points, XYZC.
%
%          NOTES:  1.  For determining the femur and tibia cartilage
%                  points that are in contact within the threshold.
%
%                  2.  The X and Y or X, Y and Z coordinates for the
%                  points must be in the columns of the coordinate
%                  arrays.
%
%                  3.  The points for the piece-wise linear line must
%                  be in order to form a piece-wise linear line.  The
%                  contacting points indexed by NID match the points
%                  on the piece-wise linear line.
%
%                  4.  M-file pts2lin.m must be in the current
%                  directory or path.

%          18-Nov-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error([' *** ERROR in PWL_CONTCT:  At least three (3) inputs', ...
         ' are required!']);
end
%
% Check Inputs
%
[nf,ncol1] = size(xyzf);
[nt,ncol2] = size(xyzt);
if ncol1~=ncol2||(ncol1~=3&&ncol1~=2)
  error([' *** ERROR in PWL_CONTCT:  All input points must have', ...
         ' two or three columns!']);
end
%
% Threshold Squared
%
thresh2 = thresh*thresh;               % Threshold squared
%
% Get Number of Line Segments and Initialize Variables
%
nl = nf-1;              % Number of line segments
%
nid = [];               % Index to points within threshold of PWL line
xyzc = [];              % Coordinates of points on PWLL within threshold
%
% Loop Along the Segments of the Piece-Wise Linear Line (PWLL)
%
for k = 1:nl
   [xyzp,t] = pts2lin(xyzf(k,:),xyzf(k+1,:)-xyzf(k,:),xyzt);
   idl = t>=0&t<=1;     % Within the line?
   if any(idl)
     d = xyzt(idl,:)-xyzp(idl,:);      % Difference
     dist = sum(d.*d,2);               % Distance
%
     idl = find(idl);
     idc = dist<=thresh2;
%
     if any(idc)        % Within threshold?
       nc = idl(idc);
       xyzc = [xyzc; xyzp(nc,:)];      % Position on PWL
       nid = [nid; nc];                % Index to points within threshold
     end
%
   end
end
%
% Check End Points Where Distance May Not Be Perpendicular
%
nid2 = [];              % Index to points within threshold of PWL line
xyzc2 = [];             % Coordinates of points on PWLL within threshold
%
for k = 1:nf
   d = xyzt-repmat(xyzf(k,:),nt,1);    % Difference
   dist = sum(d.*d,2);  % Distance
   idc = dist<=thresh2;
   nc = find(idc);
%
   if ~isempty(nc)      % Within threshold?
     idn = setdiff(nc,nid);
     if ~isempty(idn)   % Not already at minimum distance
       nid2 = [nid2; idn];
       nnew = size(nid2,1);
       xyzc2 = [xyzc2; repmat(xyzf(k,:),nnew,1)];
     end
   end
end
%
% Get Unique Minimum Points Among End Points
%
if ~isempty(nid2)
%
  [uid,inc] = unique(nid2);
  [~,inc2] = unique(nid2,'last');
  idu = inc~=inc2;      % Unique points?
%
  if any(idu)           % More than one point with same ID
%
    idnu = uid(idu);    % Get IDs with more than one point
    nnu = size(idnu);   % Number of IDs with more than one point
%
    for k = 1:nnu
%
       idx = nid2==idnu(k);
%
       d = xyzt(idx,:)-xyzc2(idx,:);   % Difference
       dist = sum(d.*d,2);             % Distance
       [~,idm] = min(dist);            % Get index to minimum
%
       idxx = find(idx);
       idm = idxx(idm);
       idx(idm) = false;
       idx = ~idx;      % Logical index to unique IDs
%
       nid2 = nid2(idx);
       xyzc2 = xyzc2(idx,:);
%
    end
  end
%
% Combine Points
%
  nid = [nid; nid2];
  [nid,ids] = sort(nid);
  xyzc = [xyzc; xyzc2];
  xyzc = xyzc(ids,:);
%
end
%
return