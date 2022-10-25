rad = 6;
%
figure;
orient landscape;
%
% Plot Cylinder
%
[xc,yc,zc] = cylinder(rad*ones(4,1),72);
zc = 45*zc-15;          % zc = 50*zc-20; for function?
plot3(xc,yc,zc,'k-');
hold on;
plot3(xc',yc',zc','k-');
%
% Loop through Slices
%
for k = 7:16
%
slk = rsl(k);
%
% Get Coordinates for Both Lines
%
fb = ft1{1}{1,k};
tb = ft1{1}{2,k};
%
% Plot Lines
%
% if iplt?
plot3(fb(:,1),fb(:,2),fb(:,3),'b.-');
hold on;
plot3(tb(:,1),tb(:,2),tb(:,3),'g.-');
axis equal;
%
plot3(fb(1,1),fb(1,2),fb(1,3),'b^');   % Start
plot3(fb(19,1),fb(19,2),fb(19,3),'bs');% End
%
plot3(tb(1,1),tb(1,2),tb(1,3),'g^');   % Start
plot3(tb(10,1),tb(10,2),tb(10,3),'gs');% End
%
% First Line (Femoral Bone) Intersections
%
nf = size(fb,1)-1;
n = zeros(nf,1);
xyz = cell(nf,1);
%
for l = 1:nf
   [n(l),xyz{l}] = cyl_lin3(fb(l,:),fb(l+1,:),rad);
   if sum(n)==2
     break;
   end
end
%
if sum(n)~=2
  continue;
end
%
xyzi = cell2mat(xyz(logical(n)));
plot3(xyzi(:,1),xyzi(:,2),xyzi(:,3),'rs');
%
idn = find(n);
if size(idn,1)==1
  idn = [];
else
  idn = idn(1)+1:idn(2);
end
%
xyzf = [xyzi(1,:); fb(idn,:); xyzi(2,:)];
%
% Second Line (Tibial Bone) Intersections
%
nt = size(tb,1)-1;
n = zeros(nt,1);
xyz = cell(nt,1);
%
for l = 1:nt
   [n(l),xyz{l}] = cyl_lin3(tb(l,:),tb(l+1,:),rad);
   if sum(n)==2
     break;
   end
end
%
if sum(n)~=2
  continue;
end
%
xyzi = cell2mat(xyz(logical(n)));
plot3(xyzi(:,1),xyzi(:,2),xyzi(:,3),'rs');
%
idn = find(n);
if size(idn,1)==1
  idn = [];
else
  idn = idn(1)+1:idn(2);
end
%
xyzt = [xyzi(1,:); tb(idn,:); xyzi(2,:)];
%
plot3(xyzt(:,1),xyzt(:,2),xyzt(:,3),'go','MarkerSize',7,'LineWidth',1);
plot3(xyzf(:,1),xyzf(:,2),xyzf(:,3),'bo','MarkerSize',7,'LineWidth',1);
%
% Check Directions of Lines - Should be Opposite Directions
%
df = xyzf(end,:)-xyzf(1,:);
dt = xyzt(end,:)-xyzt(1,:);
dchk = df*dt';          % Dot product - should be negative
if dchk>0
  xyzt = flipud(xyzt);
  dt = xyzt(end,:)-xyzt(1,:);
  dchk = df*dt';
  if dchk>=0
    error(' *** ERROR in ???: ???');
  end
end
%
% Create and Plot Polygon
%
xyzp = [xyzf; xyzt];
hp = patch(xyzp(:,1),xyzp(:,2),xyzp(:,3),'k');
set(hp,'FaceAlpha',0.5);
set(hp,'EdgeColor','none');
%
end
%
return