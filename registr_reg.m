%#######################################################################
%
%                 * REGISTRation REGression Program *
%
%          M-File which reads the registration MS-Excel spreadsheet, 
%     Results\Registration_data.xlsx, and finds extreme values of
%     registration and performs a regression between the 3D and the two
%     2D registrations.
%
%          Extreme values are written to the screen and the regression
%     plots are written to the Postscript file registr_reg.ps.
%
%     NOTES:  1.  MS-Excel spreadsheet, Results\Registration_data.xlsx,
%             must be on the specified path.
%
%             2.  M-file regres3.m must be in the current directory or
%             path.
%
%     14-Jan-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Read Registration Data MS-Excel Spreadsheet
%
xlsnam = 'Results\Registration_data.xlsx';
[~,xshts] = xlsfinfo(xlsnam);
nshts = size(xshts,2);
%
tdata = cell(nshts,1);
%
for k = 1:nshts
   t = readtable(xlsnam,'Sheet',xshts{k});
   nrow = size(t,1);
   Visit = repmat(str2num(xshts{k}(end)),nrow,1);
   Subject = repmat(xshts{k}(1:5),nrow,1);
   tdata{k} = addvars(t,Subject,Visit,'before',1);
end
%
% Get Just the Transformations
%
tdata = vertcat(tdata{:});
vnams = tdata.Properties.VariableNames;
transf = table2array(removevars(tdata,vnams(1:5)));
%
% Get Group IDs
%
tid = table2cell(removevars(tdata,vnams(6:end)));
subj = char(tid(:,1));
subjn = str2num(subj(:,1:2));
visit = logical(cell2mat(tid(:,2))-1);
type = startsWith(tid(:,4),'T2s');
%
% X and Y Translations Only
% First Two Columns are 3D
% Last Four Columns are 2D for Two Different Slices
%
xy = [transf(:,1:2) transf(:,7:10)];
x3d = xy(:,1);
y3d = xy(:,2);
x2d1 = xy(:,3);
y2d1 = xy(:,4);
x2d2 = xy(:,5);
y2d2 = xy(:,6);
%
% Find Large Translations
%
id3 = find(any((abs(xy)>3)')');
n3 = size(id3,1);
id6 = find(any((abs(xy)>6)')');
n6 = size(id6,1);
%
t3 = tdata(id3,:)
t6 = tdata(id6,:)
%
% Regression Between the 3D and the Two 2D Registrations
%
[bx1,r2_x1] = regres3(x2d1,x3d,1);     % Slice 1
[bx2,r2_x2] = regres3(x2d2,x3d,1);     % Slice 2
%
[by1,r2_y1] = regres3(y2d1,y3d,1);     % Slice 1
[by2,r2_y2] = regres3(y2d2,y3d,1);     % Slice 2
%
% Plot X Regressions
%
figure;
orient landscape;
plot(x2d1,x3d,'gx','Color',[0 0.66 0],'LineWidth',1,'MarkerSize',7);
hold on;
plot(x2d2,x3d,'bo','LineWidth',1,'MarkerSize',7);
xr = [min([x2d1; x2d2]); max([x2d1; x2d2])];
yr1 = polyval(bx1,xr);
yr2 = polyval(bx2,xr);
plot(xr,yr1,'r-','Color', [1 0.3 0.3],'LineWidth',1);
plot(xr,yr2,'r-','Color', [0.6 0 0],'LineWidth',1);
axis equal;
legend({'Slice 1','Slice 2','Best-Fit Slice 1', 'Best-Fit Slice 2'}, ...
       'Location','northeast');
xlabel('2D X Translations (pixels)','FontSize',12,'FontWeight','bold');
ylabel('3D X Translations (pixels)','FontSize',12,'FontWeight','bold');
title('Correlations between 3D and 2D Registrations','FontSize',16, ...
      'FontWeight','bold');
%
print -dpsc2 -fillpage -r600 registr_reg.ps
%
axis([-3 3 -3 3]);
strx1 = sprintf('X_{3D} = %.3f*X_{2D_{Slice 1}}%+.3f,  R^2 = %.3f', ...
                bx1,r2_x1);
strx2 = sprintf('X_{3D} = %.3f*X_{2D_{Slice 2}}%+.3f,  R^2 = %.3f', ...
                bx2,r2_x2);
text(-2.5,-1.5,strx1,'FontSize',12);
text(-2.5,-2,strx2,'FontSize',12);
%
print -dpsc2 -fillpage -r600 -append registr_reg.ps
%
% Plot Y Regressions
%
figure;
orient landscape;
plot(y2d1,y3d,'gx','Color',[0 0.66 0],'LineWidth',1,'MarkerSize',7);
hold on;
plot(y2d2,y3d,'bo','LineWidth',1,'MarkerSize',7);
xr = [min([y2d1; y2d2]); max([y2d1; y2d2])];
yr1 = polyval(by1,xr);
yr2 = polyval(by2,xr);
plot(xr,yr1,'r-','Color', [1 0.3 0.3],'LineWidth',1);
plot(xr,yr2,'r-','Color', [0.6 0 0],'LineWidth',1);
axis equal;
legend({'Slice 1','Slice 2','Best-Fit Slice 1', 'Best-Fit Slice 2'}, ...
       'Location','northeast');
xlabel('2D Y Translations (pixels)','FontSize',12,'FontWeight','bold');
ylabel('3D Y Translations (pixels)','FontSize',12,'FontWeight','bold');
title('Correlations between 3D and 2D Registrations','FontSize',16, ...
      'FontWeight','bold');
%
print -dpsc2 -fillpage -r600 -append registr_reg.ps
%
axis([-5 5 -5 5]);
stry1 = sprintf('Y_{3D} = %.3f*Y_{2D_{Slice 1}}%+.3f,  R^2 = %.3f', ...
                by1,r2_y1);
stry2 = sprintf('Y_{3D} = %.3f*Y_{2D_{Slice 2}}%+.3f,  R^2 = %.3f', ...
                by2,r2_y2);
text(-4.25,-2.5,stry1,'FontSize',12);
text(-4.25,-3.25,stry2,'FontSize',12);
%
print -dpsc2 -fillpage -r600 -append registr_reg.ps
%
% Differences between T1rho and T2*?
%
[bt0x1,r2t0_x1] = regres3(x2d1(~type),x3d(~type),1)   % T1rho Slice 1
[bt0x2,r2t0_x2] = regres3(x2d2(~type),x3d(~type),1)   % T1rho Slice 2
%
[bt1x1,r2t1_x1] = regres3(x2d1(type),x3d(type),1)     % T2* Slice 1
[bt1x2,r2t1_x2] = regres3(x2d2(type),x3d(type),1)     % T2* Slice 2
%
[bt0y1,r2t0_y1] = regres3(y2d1(~type),y3d(~type),1)   % T1rho Slice 1
[bt0y2,r2t0_y2] = regres3(y2d2(~type),y3d(~type),1)   % T1rho Slice 2
%
[bt1y1,r2t1_y1] = regres3(y2d1(type),y3d(type),1)     % T2* Slice 1
[bt1y2,r2t1_y2] = regres3(y2d2(type),y3d(type),1)     % T2* Slice 2
%
return