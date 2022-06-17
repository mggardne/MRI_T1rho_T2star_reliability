%#######################################################################
%
%                     * PLoT CSV Files 4 Program *
%
%          M-File which reads the cartilage and bone digitization CSV
%     files and plots the data for visual verification.  The program
%     assumes the femur has three sets of digitizations (lateral,
%     medial and patella), the tibia has two sets of digitizations
%     (lateral and medial) and the patella has one set of digitizations.
%     Due to the angle of the knee alignment with the scanner, the
%     patella data is ignored.
%
%     NOTES:  1.  Matlab M-files plt_datsl.m and rd_roi5.m must be in
%             the current directory or path.
%
%             2.  Plots are output to PS file:  dig_plt4.ps
%
%     09-Mar-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Directories with Digitization CSV Files
%
% fdirs = {'3.0 NSA1'; '3.0 NSA2'};
fdirs = {'T2_1MAR2021'};
%
% Loop Through Directories and Plot Digitizations
%
ido = false;
%
% for k = 1:2
for k = 1:1
%
   figure;
   orient landscape;
   fdir = fdirs{k};
%
% Patella
%
   if ido
     cdat = rd_roi5(fullfile(fdir,'PAT_CART.csv'));
     cdat = cdat.data;
     plt_datsl(cdat,'m.-');
     hold on;
     bdat = rd_roi5(fullfile(fdir,'PAT_BONE.csv'));
     bdat = bdat.data;
     plt_datsl(bdat,'k.-');
   end
%
% Femur
% 3 - TROCHLEA
%
   fcdat = rd_roi5(fullfile(fdir,'FEM_CART.csv'));
   nams = {fcdat.name}';
   fcdat1 = fcdat(1).data;
   plt_datsl(fcdat1,'r.-');
   fcdat2 = fcdat(2).data;
   plt_datsl(fcdat2,'b.-');
   if ido
     fcdat3 = fcdat(3).data;
     plt_datsl(fcdat3,'m.-');
   end
   fbdat = rd_roi5(fullfile(fdir,'FEM_BONE.csv'));
   fbdat1 = fbdat(1).data;
   fbdat2 = fbdat(2).data;
   plt_datsl(fbdat1,'k.-');
   plt_datsl(fbdat2,'k.-');
   if ido
     fbdat3 = fbdat(3).data;
     plt_datsl(fbdat3,'k.-');
   end
%
% Tibia
%
   tcdat = rd_roi5(fullfile(fdir,'TIB_CART.csv'));
   tcdat1 = tcdat(1).data;
   tcdat2 = tcdat(2).data;
   plt_datsl(tcdat1,'r.-');
   plt_datsl(tcdat2,'b.-');
   tbdat = rd_roi5(fullfile(fdir,'TIB_BONE.csv'));
   tbdat1 = tbdat(1).data;
   tbdat2 = tbdat(2).data;
   plt_datsl(tbdat1,'k.-');
   plt_datsl(tbdat2,'k.-');
%
% Finish and Print Plot
%
   view(-110,15);
   axis equal;
   hx = xlabel('\leftarrow Lateral','FontSize',12,'FontWeight','bold');
   hy = ylabel('Anterior \rightarrow','FontSize',12,'FontWeight','bold');
   st = ['Series ' fdir];
   title({st; 'Red-Lateral, Magenta-Groove and Blue-Medial Cartilage'; ...
          'Black-Bone'},'FontSize',24,'FontWeight','bold', ...
          'Interpreter','none');
   if k==1
     print -dpsc2 -r600 -fillpage dig_plt4.ps
   else
     print -dpsc2 -r600 -fillpage -append dig_plt4.ps
   end
%
end
%
return