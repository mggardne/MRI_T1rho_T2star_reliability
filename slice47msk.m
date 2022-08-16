load T1rho_S1001_rois.mat maskfrm masktrm;
msk = maskfrm(:,1,7)|maskfrm(:,2,7)|masktrm(:,1,7)|masktrm(:,2,7);
mask = reshape(msk,512,512);
mask4 = repmat(mask,1,1,4);
d3d = dat3d(mask4);
d3dt = reshape(d3d,120,4)';
vd = v(mask4);
vdt = reshape(vd,120,4)';
load img_chk
v47m = v47(mask4);
v47mt = reshape(v47m,120,4)';
d3 = d3dt-vdt;
dv = v47mt-vdt;
d = d3dt-v47mt;
figure;
orient landscape;
imagesc([d3;dv;d]);
colormap gray;
axis image;
colorbar
hold on;
hl = plot([0.5 0.5; 120.5 120.5],[4.5 8.5; 4.5 8.5],'r-','LineWidth',1);
% set(gca,'YTick',1:12);
title({'Image Differences';'1-4 DFT registered-raw image'; ...
       '5-8 elastix registered-raw image'; ...
       '9-12 DFT-elastix'},'FontSize',16,'FontWeight','bold');