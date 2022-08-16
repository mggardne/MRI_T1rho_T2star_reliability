%#######################################################################
%
%                * T1Rho Series 301 CHeck FIT Program *
%
%          M-File which reads the masks from the segmentation MAT file,
%     runs dicom_chk2.m and compares the T1Rho values with the results 
%     from mri_fitr2.m.
%
%          The comparison is for subject 1 on visit 2 for series 301.
%     Series 301 is as unloaded left leg T1rho scan.  The comparison is
%     only for slice 47 in the medial compartment.
%
%     NOTES:  1.  The default Matlab directory must be subject
%             directory "MRIR 01-BC" and visit subdirectory "Visit2".
%
%             2.  MAT file T1rho_S301_rois.mat must be in subject
%             directory "MRIR 01-BC" and visit subdirectory "Visit2".
%
%             3.  T1rho MAT file mri_fitr2.mat must be in the directory
%             "..\..\Results\NACOB_Final".
%
%             4.  The M-files T1r3d_calc.m, dftreg.m, dftregistration.m
%             and exp1_fun.m must be in the current path or directory.
%
%     29-Jun-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Get Masks for Series 301
%
load T1rho_S301_rois.mat;
%
% Get T1rho for the Medial Compartment Slices
%
dicom_chk2;             % Select Series 301
%
% Get T1rho Values from mri_fitr2.mat
%
s = load(fullfile('..\..\Results\NACOB_Final','mri_fitr2.mat'));
t1r = squeeze(s.t1r_respx(1,2,1,1,2,:,:));  
nps = squeeze(s.t1r_nps(1,2,1,1,2,:,:));
%
% Combine Femur Superficial and Deep and Tibia Superficial and Deep
% Masks
%
msk = maskfrm(:,1,7)|maskfrm(:,2,7)|masktrm(:,1,7)|masktrm(:,2,7);
%
% Get and Plot dicom_chk2 T1rho Values
%
img = T1rhonls(:,:,7);  % Slice 47 (7 of 8 medial compartment slices)
img(~msk) = NaN;        % Do not plot values not within combined mask
%
figure;
orient landscape;
%
imagesc(img,[0 70]);
%
colormap jet;
axis image;
colorbar;
axis([230 260 240 265]);
title({'T1rho for Series 301'; ['Left Leg Unloaded - Medial ', ...
       'Compartment']; 'dicom_chk values'},'FontSize',16, ...
       'FontWeight','bold','Interpreter','none');
print -dpsc2 -r600 -fillpage T1r_chk_unloaded.ps
%
% Get and Plot Slice 47 Data for the Medial Compartment
%
img2 = NaN(512);        % Blank image
idx = sum(nps{1,1}(1:6))+1:sum(nps{1,1}(1:7));   % Femur deep
img2(maskfrm(:,2,7)) = t1r{1,1}(idx);
%
idx = sum(nps{2,1}(1:6))+1:sum(nps{2,1}(1:7));   % Tibia deep
img2(masktrm(:,2,7)) = t1r{2,1}(idx);
%
idx = sum(nps{1,2}(1:6))+1:sum(nps{1,2}(1:7));   % Femur superficial
img2(maskfrm(:,1,7)) = t1r{1,2}(idx);
%
idx = sum(nps{2,2}(1:6))+1:sum(nps{2,2}(1:7));   % Tibia superficial
img2(masktrm(:,1,7)) = t1r{2,2}(idx);
%
figure;
orient landscape;
%
imagesc(img2,[0 70]);
%
colormap jet;
axis image;
colorbar;
axis([230 260 240 265]);
title({'T1rho for Series 301'; ['Left Leg Unloaded - Medial ', ...
       'Compartment']; 'mri_fitr2 values'},'FontSize',16, ...
       'FontWeight','bold','Interpreter','none');
print -dpsc2 -r600 -fillpage -append T1r_chk_unloaded.ps
%
% Get T1rho Values
%
% Column 1 - Reference T1rho values from dicom_chk2.m
% Column 2 - T1rho values from mri_fitr2.m
% Column 3 - Differences in T1rho values (dicom_chk2.m - mri_fitr2.m)
%
idv = ~isnan(img2);
t1rv = [img(idv) img2(idv) img(idv)-img2(idv)];
%
% Histogram of Differences
%
figure;
orient landscape;
%
histogram(t1rv(:,3));
%
title({'Histogram of T1\rho Differences for Series 301'; ...
       'Left Leg Unloaded - Medial Compartment'},'FontSize',16, ...
       'FontWeight','bold');
xlabel('T1\rho Differences (ms)','FontSize',12,'FontWeight','bold');
ylabel('Frequency','FontSize',12,'FontWeight','bold');
print -dpsc2 -r600 -fillpage -append T1r_chk_unloaded.ps
%
% Get Masks for All Eight (8) Slices in the Medial Compartment
%
mska = false(npx*npx,8);
t1rr = cell(8,1);
%
for k = 1:8
   mska(:,k) = maskfrm(:,1,k)|maskfrm(:,2,k)|masktrm(:,1,k)|masktrm(:,2,k);
   t1rrt = T1rhonls(:,:,k);
   t1rr{k} = t1rrt(mska(:,k));
end
%
% Plot T1rho Values
%
t1rr = cell2mat(t1rr);                 % All eight slices - dicom_chk2.m
t1ra = cell2mat(t1r(:));               % All eight slices - mri_fitr2.m
%
figure;
orient landscape;
%
plot(t1rr);             % dicom_chk2.m values
hold on;
plot(t1ra);             % mri_fitr2.m values
%
title({'T1\rho for Series 301'; ['Left Leg Unloaded - Medial ', ...
       'Compartment']},'FontSize',16,'FontWeight','bold');
ylabel('T1\rho (ms)','FontSize',12,'FontWeight','bold');
xlabel('Pixel','FontSize',12,'FontWeight','bold');
legend({'dicom_chk','mri_fitr2'},'Interpreter','none');
print -dpsc2 -r600 -fillpage -append T1r_chk_unloaded.ps
%
return