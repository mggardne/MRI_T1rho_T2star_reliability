function [tri,xy,nt] = mk_tri_2d(dat,tol,iplt)
%MK_TRI_2D Makes a triangular mesh using boundary line data from
%        a two-dimensional digitized MRI slice.  The first line is
%        assumed to be cartilage and the second line is assumed to be
%        bone.
%
%        [TRI,XY] = MK_TRI_2D(DAT) given a cell array containing two (2)
%        columns matrices with boundary line coordinate point data,
%        DAT, returns a three (3) column triangle connectivity matrix,
%        TRI and X and Y coordinates in a two columns matrix XY.
%
%        [TRI,XY,NT] = MK_TRI_2D(DAT) returns the number of triangles,
%        NT.
%
%        NOTES:  1.  Each boundary coordinate data matrix must
%                correspond to one index into the cell array DAT.
%
%                2.  The coordinates should be ordered in the same
%                direction in each line.  The dot product of the
%                directions of the adjacent lines are used to check the
%                ordering direction and the ordering direction is
%                reversed if the dot product is negative.
%
%                3.  The second boundary line may be truncated if it
%                extends pass the first boundary line.
%
%                4.  M-files lsect2.m and lsect2a.m must be in the
%                current directory or path.
%
%                5.  Similar to mk_tri_2df.m except the slopes of the
%                end points are in the opposite directions.
%
%        04-Nov-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<1)
  error(' *** ERROR in mk_tri_2d:  No input data!');
end
%
if nargin<2
  tol = 0.1;
end
if isempty(tol)
  tol = 0.1;
end
%
if nargin<3
  iplt = false;
end
%
dat = dat(:);
nslice = size(dat,1);
%
if nslice~=2
  error([' *** ERROR in mk_tri_2d:  Input cell array must have ', ...
         'two elements containing 2D coordinates for two lines!']);
end
%
% Get First (Cartilage) Line and Slopes at the Ends of the Line
%
xy1 = dat{1};
npts1 = size(xy1,1);
%
vec1 = xy1(npts1,:)-xy1(1,:);          % Direction of line
vec1 = vec1./norm(vec1);
%
dd = diff(xy1);
de = dd([1,npts1-1],:);                % Slopes at ends of top line
%
% Ends of First (Cartilage) Line
%
mp = -de(:,1)./de(:,2); % 90 degrees
ids = abs(de(:,2))<1e-8;
mp(ids) = sign(mp(ids))*1e+4;          % Large but not infinite slope
% mp(2,1) = (de(2,1)+de(2,2))./(de(2,1)-de(2,2));  % 45 degrees
% mp(1,1) = (de(1,2)-de(1,1))./(de(1,1)+de(1,2));  % 45 degrees
% ratio = tan(pi/3);      % Opposite to adjacent for 60 degrees (pi/3 radians)
% mp(2,1) = (de(2,2)+ratio*de(2,1))./(de(2,1)-ratio*de(2,2));     % 60 degrees
% mp(1,1) = (ratio*de(1,2)-de(1,1))./(de(1,1)+ratio*de(1,2));     % 60 degrees
xp = xy1([1; npts1],1);
yp = xy1([1; npts1],2);
bp = yp-mp.*xp;
%
% Get Second (Bone) Line
%
xy2 = dat{2};
npts2 = size(xy2,1);
%
vec2 = xy2(npts2,:)-xy2(1,:);
vec2 = vec2./norm(vec2);
%
% Check for Slices with a Reverse Digitization
%
dotp = vec1*vec2';
%
if dotp<tol
  xy2 = flipud(xy2);
  vec = xy2(npts2,:)-xy2(1,:);
  vec = vec./norm(vec);
  dotp2 = vec1*vec';
  if dotp2<dotp         % Revert back to original ordering
    warning([' *** WARNING in mk_tri_2d:  Ordering of points', ...
             ' in the slices may not be in the same direction!']);
    xy2 = flipud(xy2);
  end
end
%
% Cut Off Extra Points on Second (Bone) Line
% Fails if Second (Bone) is Not Longer Than First (Cartilage) Line
%
[~,id1] = lsect2a(mp(1),bp(1),xy2);
if size(id1,1)>1
  id1 = max(id1);
end
chk = max([ceil(npts2/1.7); (npts2-npts1)+1]);
if isempty(id1)||id1>chk
  id1 = 1;
end
%
[~,id2] = lsect2a(mp(2),bp(2),xy2);
if size(id2,1)>1
  id2 = min(id2);
end
if isempty(id2)||id2<npts2-chk
  id2 = npts2-1;
end
%
idc = id1:id2+1;
nptc = length(idc);
xy2 = xy2(idc,:);
%
xy = [xy1; xy2];
%
% Delaunay Triangulation
%
n = [0; cumsum([npts1; nptc])];
%
c1 = [(n(1)+1:n(2)-1)' (n(1)+2:n(2))'; n(2) n(3); (n(3):-1:n(2)+2)' ...
      (n(3)-1:-1:n(2)+1)'; n(2)+1 n(1)+1];       % Constraints
%
dt1 = delaunayTriangulation(xy,c1);
idin = isInterior(dt1);
tri = dt1(idin,:);
nt = size(tri,1);
%
% Plot Triangulations?
%
if iplt
%
  h1 = figure;
  orient tall;
%
  xt = xy(:,1);
  yt = xy(:,2);
%
  plot(xt,yt,'k.','LineWidth',1,'MarkerSize',7);
  hold on;
  npts = size(xt,1);
  text(xt,yt,int2str((1:npts)'),'Color','k','FontSize',10);
%
  trimesh(tri,xt,yt);
  text(mean(xt(tri),2),mean(yt(tri),2),int2str((1:nt)'), ...
       'Color','r','FontSize',10);
%
  h2 = figure;
  orient tall;
  plot(xy2(:,1),xy2(:,2),'k.-','LineWidth',1,'MarkerSize',7);
  hold on;
  plot(xy1(:,1),xy1(:,2),'b.-','LineWidth',1,'MarkerSize',7);
  text(xt,yt,int2str((1:npts)'),'Color','k', ...
       'FontSize',10);
%
  xp = reshape(xy(tri,1),nt,3)';
  yp = reshape(xy(tri,2),nt,3)';
  xp = repmat(mean(xp),3,1)+0.75*(xp-repmat(mean(xp),3,1));
  yp = repmat(mean(yp),3,1)+0.75*(yp-repmat(mean(yp),3,1));
  patch(xp,yp,[1 0.7 0.7]);
  text(mean(xt(tri),2),mean(yt(tri),2),int2str((1:nt)'), ...
       'Color','r','FontSize',10);
%
  axis equal;
  pause;
  close(h1,h2);
%
end
%
return