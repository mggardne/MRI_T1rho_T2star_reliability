%
fdat = [ -0.2500    2.7500
          0.7500    1.7500
          1.7500    0.7500
          2.7500   -0.2500
          3.7500    0.7500
          4.7500    1.7500
          5.7500    2.7500 ];
%
tdat = [  0.2500    0.2500
          1.2500    1.2500
          2.2500    2.2500
          3.2500    3.2500
          4.2500    2.2500
          5.2500    1.2500
          6.2500    0.2500 ];
%
[dis,idd,xyzd] = dis2lin(fdat,tdat);
%
figure;
plot(fdat(:,1),fdat(:,2),'b.-','MarkerSize',15);
hold on;
axis equal;
plot(tdat(:,1),tdat(:,2),'g.-','MarkerSize',15);
plot([xyzd(idd,1) tdat(:,1)]',[xyzd(idd,2) tdat(:,2)]','r-');
plot(xyzd(idd,1),xyzd(idd,2),'rs','MarkerSize',12);
%