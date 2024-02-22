function h = ax_arrow2(vlen,r,axdir,offst,alen)
%AX_ARROW2 Plots cylinders as axis arrows for three directions.
%
%         H = AX_ARROW2(VLEN,R,AXDIR) returns the handles H to the
%         cylinders given the arrow lengths, VLEN, arrow stem radius, R
%         and arrow directions in a three (3) by three (3) matrix.  If
%         VLEN is a scalar, all the arrows have the same length.  AXDIR
%         contains the three-dimensional directions of the arrows in
%         the rows of the matrix.  H is a three by two matrix of
%         handles with the first column is the stem of the three arrows
%         and the second column is the head of the arrow.
%
%         AX_ARROW2(VLEN,R,AXDIR,OFFST,ALEN) offsets the three arrows
%         from the origin to the three-dimensional coordinate, OFFST.
%         The length of the arrow heads is in ALEN.  Default length is
%         seven (7) times the arrow stem radius.
%
%         27-Apr-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<3
  error(' *** ERROR:  AX_ARROW2 requires three input arguments.');
end
%
if size(vlen(:),1)<3
  vlen = repmat(vlen(1),3,1);
end
%
if nargin<4||isempty(offst)
  offst = zeros(3,1);
end
%
offst = offst(:);
%
if nargin<5
  alen = 7*r;
end
%
% Make Direction Vectors Unit Vectors
%
for k = 1:3
   axdir(k,:) = axdir(k,:)/norm(axdir(k,:));
end
%
% Translation of Arrow Head
%
th = vlen(:)-alen(:);
%
% Create Arrow Stem and Head Cylinders
%
[x3,y3,z3] = cylinder(repmat(r,4,1),36);         % Stem
[xh3,yh3,zh3] = cylinder(2*r:-2*r/4:0,36);       % Head
xh3 = [zeros(1,37); xh3];
yh3 = [zeros(1,37); yh3];
zh3 = [zh3(1,:); zh3];
zh3 = zh3*alen;
%
x1 = (vlen(1)-0.99*alen)*z3;
y1 = x3;
z1 = y3;
%
xh1 = zh3+th(1);
yh1 = xh3;
zh1 = yh3;
%
x2 = x3;
y2 = (vlen(2)-0.99*alen)*z3;
z2 = y3;
%
xh2 = xh3;
yh2 = zh3+th(2);
zh2 = yh3;
%
z3 = (vlen(3)-0.99*alen)*z3;
zh3 = zh3+th(3);
%
% Rotate and Translate Arrows
%
dim = size(x1);        % Dimensions of arrow stem coordinate matrices
%
a1 = [x1(:)'; y1(:)'; z1(:)'];
a1 = axdir'*a1;
a1 = a1+repmat(offst,1,prod(dim));
x1 = reshape(a1(1,:),dim);
y1 = reshape(a1(2,:),dim);
z1 = reshape(a1(3,:),dim);
%
a2 = [x2(:)'; y2(:)'; z2(:)'];
a2 = axdir'*a2;
a2 = a2+repmat(offst,1,prod(dim));
x2 = reshape(a2(1,:),dim);
y2 = reshape(a2(2,:),dim);
z2 = reshape(a2(3,:),dim);
%
a3 = [x3(:)'; y3(:)'; z3(:)'];
a3 = axdir'*a3;
a3 = a3+repmat(offst,1,prod(dim));
x3 = reshape(a3(1,:),dim);
y3 = reshape(a3(2,:),dim);
z3 = reshape(a3(3,:),dim);
%
dim = size(xh1);        % Dimensions of arrow hat coordinate matrices
%
a1 = [xh1(:)'; yh1(:)'; zh1(:)'];
a1 = axdir'*a1;
a1 = a1+repmat(offst,1,prod(dim));
xh1 = reshape(a1(1,:),dim);
yh1 = reshape(a1(2,:),dim);
zh1 = reshape(a1(3,:),dim);
%
a2 = [xh2(:)'; yh2(:)'; zh2(:)'];
a2 = axdir'*a2;
a2 = a2+repmat(offst,1,prod(dim));
xh2 = reshape(a2(1,:),dim);
yh2 = reshape(a2(2,:),dim);
zh2 = reshape(a2(3,:),dim);
%
a3 = [xh3(:)'; yh3(:)'; zh3(:)'];
a3 = axdir'*a3;
a3 = a3+repmat(offst,1,prod(dim));
xh3 = reshape(a3(1,:),dim);
yh3 = reshape(a3(2,:),dim);
zh3 = reshape(a3(3,:),dim);
%
% Plot Arrows
%
h(3,2) = surfl(xh3,yh3,zh3);
hold on;
h(3,1) = surfl(x3,y3,z3);
%
h(2,2) = surfl(xh2,yh2,zh2);
h(2,1) = surfl(x2,y2,z2);
%
h(1,2) = surfl(xh1,yh1,zh1);
h(1,1) = surfl(x1,y1,z1);
%
set(h,'EdgeColor','none','FaceColor','r');
%
return