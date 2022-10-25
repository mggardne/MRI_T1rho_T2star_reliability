r = 6;
npts = 24;
t = (0:pi/36:2*pi)';
x = r*cos(t);
y = r*sin(t);
%
ids = [4 32 32 68 42 68;
       32 4 68 32 68 42];
%
tsa = false(2,6);
dts = zeros(2,6);
%
for k = 1:6
%
   figure;
   plot(x,y,'r.-');
   axis equal;
   hold on;
   text(x,y,int2str((1:73)'),'Color','r');
%
   endpts = [x(ids(:,k)) y(ids(:,k))];
   xe = endpts(:,1);
   ye = endpts(:,2);
   te = acos(xe./r);
   ts = asin(ye./r);
   ts = sign(ts)<0;
   tsa(:,k) = ts;
   te(ts) = 2*pi-te(ts);     % Extend range to 2*pi (360 degrees)
%
   dt1 = (te(2)-te(1))/(npts+1);
   dts(1,k) = sign(dt1);
   t1 = (te(1):dt1:te(2))';
   x1 = r*cos(t1);
   y1 = r*sin(t1);
%
   h1 = plot(x1,y1,'bo-');
%
   dt2 = (te(2)-sign(dt1)*2*pi-te(1))/(npts+1);
   dts(2,k) = sign(dt2);
   t2 = (te(1):dt2:(te(2)-sign(dt1)*2*pi))';
   x2 = r*cos(t2);
   y2 = r*sin(t2);
%
   h2 = plot(x2,y2,'go-');
%
end
