%#######################################################################
%
%           * Tibia Cartilage THicKness Comparison Program *
%
%          M-File which reads the T1FFE and T1rho cartilage thicknesses,
%     finds a combined grid that covers both data sets and finds the
%     differences in cartilage thicknesses.
%
%          Plots, Outputs, ...
%
%     NOTES:  1.  Both grids must have integer coordinates.
%
%     17-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Get Analysis Grids
%
ffeg = load(fullfile('TibiaCartThk_1_14Mar2022(FFE)','tgrid08_1.mat'));
rhog = load(fullfile('TibiaCartThk_1_11Mar2022(T1Rho)', ...
            'tgrid08_1.mat'));
%
% Get Coordinates and Ranges for Lateral Compartment
%
xqlf = ffeg.xql;
yqlf = ffeg.yql;
%
xqlr = rhog.xql;
yqlr = rhog.yql;
%
xmnf = min(xqlf);
xmxf = max(xqlf);
ncf = ffeg.nxl;         % Number of columns
xmnr = min(xqlr);
xmxr = max(xqlr);
ncr = rhog.nxl;         % Number of columns
%
ymnf = min(yqlf);
ymxf = max(yqlf);
nrf = ffeg.nyl;         % Number of rows
ymnr = min(yqlr);
ymxr = max(yqlr);
nrr = rhog.nyl;         % Number of rows
%
% Get Combined Grid for Both Data Sets
%
if xmnf<xmnr
  xmn = xmnf;
else
  xmn = xmnr;
end
if xmxf>xmxr
  xmx = xmxf;
else
  xmx = xmxr;
end
%
if ymnf<ymnr
  ymn = ymnf;
else
  ymn = ymnr;
end
if ymxf>ymxr
  ymx = ymxf;
else
  ymx = ymxr;
end
%
[xgl,ygl] = meshgrid(xmn:xmx,ymn:ymx);
[nlr,nlc] = size(xgl);
xgl = xgl(:);
ygl = ygl(:);
quadl = quadconn(nlr,nlc);   % Quadrilateral connectivity for grid
%
% Get Indexes into Combined Grid
%
nl = nlr*nlc;           % Number of points in combined lateral grid
idx = reshape(1:nl,nlr,nlc);
%
offstcf = round(xmnf-xmn+1);      % Differences in integer grid == column index
offstrf = round(ymnf-ymn+1);      % Differences in integer grid == row index
idxf = idx(offstrf:offstrf+nrf-1,offstcf:offstcf+ncf-1);   % Index for T1FFE
%
offstcr = round(xmnr-xmn+1);      % Differences in integer grid == column index
offstrr = round(ymnr-ymn+1);      % Differences in integer grid == row index
idxr = idx(offstrr:offstrr+nrr-1,offstcr:offstcr+ncr-1);   % Index for T1rho
%
% Read Cartilage Thicknesses
%
ffe = load(fullfile('TibiaCartThk_1_14Mar2022(FFE)','tcart08_1.mat'));
rho = load(fullfile('TibiaCartThk_1_11Mar2022(T1Rho)','tcart08_1.mat'));
%
cthkfl = ffe.cthkls;
cthkrl = rho.cthkls;
%
% Get Cartilage Thicknesses at Similar Coordinates
%
cthklf = NaN(nlr,nlc);  % NaN == missing data
cthklr = NaN(nlr,nlc);
%
cthklf(idxf) = cthkfl;
cthklr(idxr) = cthkrl;
%
% Get Coordinates and Ranges for Medial Compartment
%
xqmf = ffeg.xqm;
yqmf = ffeg.yqm;
%
xqmr = rhog.xqm;
yqmr = rhog.yqm;
%
xmnf = min(xqmf);
xmxf = max(xqmf);
ncf = ffeg.nxm;         % Number of columns
xmnr = min(xqmr);
xmxr = max(xqmr);
ncr = rhog.nxm;         % Number of columns
%
ymnf = min(yqmf);
ymxf = max(yqmf);
nrf = ffeg.nym;         % Number of rows
ymnr = min(yqmr);
ymxr = max(yqmr);
nrr = rhog.nym;         % Number of rows
%
% Get Combined Grid for Both Data Sets
%
if xmnf<xmnr
  xmn = xmnf;
else
  xmn = xmnr;
end
if xmxf>xmxr
  xmx = xmxf;
else
  xmx = xmxr;
end
%
if ymnf<ymnr
  ymn = ymnf;
else
  ymn = ymnr;
end
if ymxf>ymxr
  ymx = ymxf;
else
  ymx = ymxr;
end
%
[xgm,ygm] = meshgrid(xmn:xmx,ymn:ymx);
[nmr,nmc] = size(xgm);
xgm = xgm(:);
ygm = ygm(:);
quadm = quadconn(nmr,nmc);   % Quadrilateral connectivity for grid
%
% Get Indexes into Combined Grid
%
nm = nmr*nmc;           % Number of points in combined medial grid
idx = reshape(1:nm,nmr,nmc);
%
offstcf = round(xmnf-xmn+1);      % Differences in integer grid == column index
offstrf = round(ymnf-ymn+1);      % Differences in integer grid == row index
idxf = idx(offstrf:offstrf+nrf-1,offstcf:offstcf+ncf-1);   % Index for T1FFE
%
offstcr = round(xmnr-xmn+1);      % Differences in integer grid == column index
offstrr = round(ymnr-ymn+1);      % Differences in integer grid == row index
idxr = idx(offstrr:offstrr+nrr-1,offstcr:offstcr+ncr-1);   % Index for T1rho
%
% Get Cartilage Thicknesses
%
cthkfm = ffe.cthkms;
cthkrm = rho.cthkms;
%
% Get Cartilage Thicknesses at Similar Coordinates
%
cthkmf = NaN(nmr,nmc);  % NaN == missing data
cthkmr = NaN(nmr,nmc);
%
cthkmf(idxf) = cthkfm;
cthkmr(idxr) = cthkrm;
%
% Calculate Differences
%
cthkld = cthklr-cthklf;
cthkmd = cthkmr-cthkmf;
idvl = ~isnan(cthkld);  % Valid differences
idvm = ~isnan(cthkmd);  % Valid differences
%
% Plots of Cartilage Thicknesses and Thickness Differences
%
figure;
orient tall;
colormap jet;
%
subplot(3,1,1);
cthk = cthklr;
cthk(~idvl) = NaN;
patch(xgl(quadl'),ygl(quadl'),cthk(quadl'),'FaceColor','interp', ...
      'EdgeColor','interp');
hold on;
cthk = cthkmr;
cthk(~idvm) = NaN;
patch(xgm(quadm'),ygm(quadm'),cthk(quadm'),'FaceColor','interp', ...
      'EdgeColor','interp');
view(-90,90);
title('T1\rho Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
colorbar;
%
subplot(3,1,2);
cthk = cthklf;
cthk(~idvl) = NaN;
patch(xgl(quadl'),ygl(quadl'),cthk(quadl'),'FaceColor','interp', ...
      'EdgeColor','interp');
hold on;
cthk = cthkmf;
cthk(~idvm) = NaN;
patch(xgm(quadm'),ygm(quadm'),cthk(quadm'),'FaceColor','interp', ...
      'EdgeColor','interp');
view(-90,90);
title('T1FFE Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
colorbar;
%
subplot(3,1,3);
dthkl = cthkld;
patch(xgl(quadl'),ygl(quadl'),dthkl(quadl'),'FaceColor','interp', ...
      'EdgeColor','interp');
hold on;
dthkm = cthkmd;
patch(xgm(quadm'),ygm(quadm'),dthkm(quadm'),'FaceColor','interp', ...
      'EdgeColor','interp');
view(-90,90);
title('T1\rho - T1FFE Cartilage Thickness Differences', ...
      'FontSize',16,'FontWeight','bold');
colorbar;
%
print -dpsc2 -r600 -fillpage tcthk_cmp.ps;
%
% Statistics and Plot of Differences
%
dthkl = dthkl(idvl(:));
dthkm = dthkm(idvm(:));
dmeanl = mean(dthkl);
dmeanm = mean(dthkm);
%
dstdl = std(dthkl);
dstdm = std(dthkm);
%
figure;
orient landscape;
%
subplot(2,1,1);
plot(dthkl,'k.');
hold on;
axlim = axis;
plot(axlim(1:2),[0 0],'k-');
plot(axlim(1:2),[dmeanl dmeanl],'k-','LineWidth',1);
plot(axlim(1:2),[dmeanl+2*dstdl dmeanl+2*dstdl],'k:','LineWidth',1);
plot(axlim(1:2),[dmeanl-2*dstdl dmeanl-2*dstdl],'k:','LineWidth',1);
plot(axlim(1:2),[dmeanl+3*dstdl dmeanl+3*dstdl],'r--','LineWidth',1);
plot(axlim(1:2),[dmeanl-3*dstdl dmeanl-3*dstdl],'r--','LineWidth',1);
title('Lateral Compartment Cartilage Thickness Differences', ...
      'FontSize',16,'FontWeight','bold');
%
subplot(2,1,2);
plot(dthkm,'k.');
hold on;
axlim = axis;
plot(axlim(1:2),[0 0],'k-');
plot(axlim(1:2),[dmeanm dmeanm],'k-','LineWidth',1);
plot(axlim(1:2),[dmeanm+2*dstdm dmeanm+2*dstdm],'k:','LineWidth',1);
plot(axlim(1:2),[dmeanm-2*dstdm dmeanm-2*dstdm],'k:','LineWidth',1);
plot(axlim(1:2),[dmeanm+3*dstdm dmeanm+3*dstdm],'r--','LineWidth',1);
plot(axlim(1:2),[dmeanm-3*dstdm dmeanm-3*dstdm],'r--','LineWidth',1);
title('Medial Compartment Cartilage Thickness Differences', ...
      'FontSize',16,'FontWeight','bold');
%
print -dpsc2 -r600 -fillpage -append tcthk_cmp.ps;
%
return