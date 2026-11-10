%#######################################################################
%
%            * MRI Reliability CARTilage SURFace Program *
%
%          M-File which reads the MRIR PTOA ROI MAT files (T*_prois.mat)
%     and gets the femur slice segmentation region designations
%     (lateral, medial, and trochlea) for the bone (idregb) and
%     cartilage (idregc) points.  These two variables, the rotation
%     matrices and translation vectors between the MRI and bone
%     coordinate systems for the femur and tibia, and the points and
%     normal vectors for the four dividing planes defining the six (6)
%     regions of interest are saved in MAT files.
%
%     NOTES:  1.  The registered ROI MAT files must be in subdirectories
%             "Visit1" or "Visit2" in subject directories "MRIR *"
%             where "*" is the subject number.
%
%             2.  For T1rho unloaded trials on visits 1 and 2.
%
%             3.  T1rho PTOA ROI MAT files must start with "T1rho_S"
%             and end with "01_prois.mat".
%
%             4.  Must start in the MRI_Reliability_Study directory.
%
%     06-Oct-2026 * Mack Gardner-Morse
%

%#######################################################################
%
% Output Directory
%
resdir = fullfile('Results','HSS');    % Results directory
%
% Get Subject Directories
%
sdirs = { 'MRIR 01-BC'
          'MRIR 02-ES'
          'MRIR 05-SC'
          'MRIR 06-MM'
          'MRIR 07-BP'
          'MRIR 08-AS'
          'MRIR 09-CR'
          'MRIR 10-JF'
          'MRIR 12-AO'
          'MRIR 13-AK' };
nsubj = size(sdirs,1);
%
% Visits
%
vdirs = {'Visit1';'Visit2'};
vds = {'V1';'V2'};
nvisit = size(vdirs,1);
%
% T1rho ROI File Names
%
rdir = 'T1rho_S*01_prois.mat';
%
% Unload Trial Only
%
ld = 0;
ldstr = 'Unloaded';
%
% Loop Through Subjects
%
for ks = 1:nsubj
%
   sdir = sdirs{ks};
%
% Loop Through Visits
%
   for kv = 1:nvisit    % Both visits
%
      vdir = vdirs{kv};                % Get visit
      vd = vds{kv};                    % Short version of visit
      pdir = fullfile(sdir,vdir);      % Get path to subject and visit
      fprintf(1,'Working on %s (%i of 10) ...\n',pdir,ks);
%
% Get ROI MAT Files in Subject at this Visit
%
      roinams = dir(fullfile(pdir,rdir));
      roinams = {roinams.name}';
      nrois = size(roinams,1);
%
% Loop Through ROI Files
%
      for km = 1:nrois
%
% Read ROI File
%
         roinam = roinams{km};
         load(fullfile(pdir,roinam),'brois','f3f','leg','plan_nvs', ...
              'plan_pts','rmat','rslf','rtmat','subj', ...
              'subjs','tplan_nvs','tplan_pts','vnum','vstr', ...
              'xyzc','xyzt');
%
         resnam = fullfile(resdir,[sdir '_' vd '_T1R_' leg]);
%
% Femur Cartilage Segmentations
%
% 1 - Lateral
% 2 - Medial
% 3 - Trochlea
%
         idlat = brois(1).rois(1).roi(1).imageno;
         idmed = brois(1).rois(1).roi(2).imageno;
         if size(brois(1).rois(1).roi,2)>2
           idtroch = brois(1).rois(1).roi(3).imageno;
         else
           idtroch = 0;
         end
         idl = ismember(rslf,idlat);
         idm = ismember(rslf,idmed);
         idt = ismember(rslf,idtroch);
%
% Femur Cartilage Segmentations Points
%
         f3fc = f3f(1,:)';
         npts = cellfun('length',f3fc);
         idregc = zeros(sum(npts),1);
%
         if idl(1)
           id1 = sum(npts(idl));
           idregc(1:id1) = 1;
           id2 = sum(npts(idt));
           idregc(id1+1:id1+id2) = 3;
           id3 = sum(npts(idm));
           idregc(id1+id2+1:id1+id2+id3) = 2;
         else
           id1 = sum(npts(idm));
           idregc(1:id1) = 2;
           id2 = sum(npts(idt));
           idregc(id1+1:id1+id2) = 3;
           id3 = sum(npts(idl));
           idregc(id1+id2+1:id1+id2+id3) = 1;
         end
%
% Femur Bone Segmentations
%
% 1 - Lateral
% 2 - Medial
% 3 - Trochlea
%
         idlat = brois(1).rois(2).roi(1).imageno;
         idmed = brois(1).rois(2).roi(2).imageno;
         if size(brois(1).rois(2).roi,2)>2
           idtroch = brois(1).rois(2).roi(3).imageno;
         else
           idtroch = 0;
         end
         idl = ismember(rslf,idlat);
         idm = ismember(rslf,idmed);
         idt = ismember(rslf,idtroch);
%
% Femur Bone Segmentations Points
%
         f3fb = f3f(2,:)';
         npts = cellfun('length',f3fb);
         idregb = zeros(sum(npts),1);
%
         if idl(1)
           id1 = sum(npts(idl));
           idregb(1:id1) = 1;
           id2 = sum(npts(idt));
           idregb(id1+1:id1+id2) = 3;
           id3 = sum(npts(idm));
           idregb(id1+id2+1:id1+id2+id3) = 2;
         else
           id1 = sum(npts(idm));
           idregb(1:id1) = 2;
           id2 = sum(npts(idt));
           idregb(id1+1:id1+id2) = 3;
           id3 = sum(npts(idl));
           idregb(id1+id2+1:id1+id2+id3) = 1;
         end
%
% Save MAT File
%
         save([resnam '.mat'],'idregb','idregc','ld','ldstr','leg', ...
              'plan_nvs','plan_pts','rmat','rtmat','subj','subjs', ...
              'tplan_nvs','tplan_pts','vnum','vstr','xyzc','xyzt');
%
      end
%
   end
%
end