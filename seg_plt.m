sdir = 'MRIR 07-BP';
vdir = 'Visit2';
%
% Tibia
%
rdir = fullfile(sdir,vdir,'RHO','Tibia');
%
fnames = {'007_R_SAG_TIB_RHO_LD_V2_CF.csv'
          '007_R_SAGAR_TIB_RHO_UL_V2_CF.csv'
          '007_R_SAGAR_TIB_RHO_LD_V2_CF.csv'};
%
fnams = strrep(fnames,'.csv','_org.csv');
%
for k = 1:3
   roi = rd_roi6(fullfile(rdir,fnams{k}));
   datl = roi(1).data;
   datm = roi(2).data;
   figure;
   plt_datsl(datl,'k.-');
   hold on;
   plt_datsl(datm,'b.-');
   roi = rd_roi6(fullfile(rdir,fnames{k}));
   datl = roi(1).data;
   datm = roi(2).data;
   plt_datsl(datl,'g.-');
   hold on;
   plt_datsl(datm,'r.-');
   axis equal;
end
%
% Femur
%
rdir = fullfile(sdir,vdir,'RHO','Femur');
%
fnames = {'007_R_SAGAR_FEM_RHO_LD_V2_JF.csv'
          '007_R_SAGAR_FEM_RHO_UL_V2_JF.csv'};
%
fnams = strrep(fnames,'.csv','_org.csv');
%
for k = 1:2
   roi = rd_roi6(fullfile(rdir,fnams{k}));
   datl = roi(1).data;
   datm = roi(2).data;
   figure;
   plt_datsl(datl,'k.-');
   hold on;
   plt_datsl(datm,'b.-');
   roi = rd_roi6(fullfile(rdir,fnames{k}));
   datl = roi(1).data;
   datm = roi(2).data;
   plt_datsl(datl,'g.-');
   hold on;
   plt_datsl(datm,'r.-');
   axis equal;
end
%
return