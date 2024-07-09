%#######################################################################
%
%                     * MRI TREAD CHecK Program *
%
%          M-file which reads the registered MRI data, segmentation,
%     and meniscus MAT files and checks to be sure the segmentation
%     slices include all the meniscus slices.
%
%     NOTES:  1.  Data MAT files must be in subject directories starting
%             with "MRIR" and visit subdirectories "Visit1" or "Visit2".
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "2rois2".  See rd_dicom.m and
%             seg_2rois2.m.  Meniscus MAT files must contain
%             "_mrois.mat".
%
%             3.  M-file to read the meniscus segmentation MAT files
%             for the slices to use for the tread analysis. 
%
%     03-Jul-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Get Subject Directories and Visit Subdirectories
%
sdirs = dir('MRIR *');
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
nsubj = size(sdirs,1);
%
vdirs = {'Visit1'; 'Visit2'};          % Visit directories
nvisit = size(vdirs,1);
%
fprintf(1,'\n\n');
%
% Loop through Subjects
%
for ks = 1:nsubj
%
% Get Subject Directory, Name and Number
%
   sdir = sdirs{ks};                   % Current subject directory
   subjnam = sdir(6:end);              % Subject name as text
   subj = eval(subjnam(1:2));          % Subject number
%
% Loop through Visits
%
   for kv = 1:nvisit
%
% Get Visit Subdirectory, Name and Number
%
      vdir = vdirs{kv};                % Current visit directory
      vstr = int2str(kv);              % Visit number as text
      vnam = ['Visit ' vstr];          % Visit name as text
      vid = kv-1;                      % Visit number
%
% Directory with Data Matlab MAT Files
%
      rdir = fullfile(sdir,vdir);      % Directory with data
%
% Get T1rho MAT Files in Directory
%
      d = dir(fullfile(rdir,'T1rho_S*.mat'));
      roinams = {d.name}';
      idr = contains(roinams,'roi','IgnoreCase',true);     % Masks
%
      rhonams = roinams(~idr);         % Image MAT files
      idv = ~contains(rhonams,'_3');   % Ignore test MAT files
      rhonams = rhonams(idv);
      nrho = size(rhonams,1);
%
      idm = contains(roinams,'_mrois.mat','IgnoreCase',true);   % Meniscus ROIs
      mroinams = roinams(idm);         % ROI MAT files
      nmroi = size(mroinams,1);
%
      idr = contains(roinams,'2rois2.mat','IgnoreCase',true);   % Small ROIs
      roinams = roinams(idr);          % ROI MAT files
      nroi = size(roinams,1);
%
      if (nrho~=nroi)||(nrho~=nmroi)
        error([' *** ERROR in mri_tread_chk:  Number of T1rho MAT', ...
               ' files does not match the number of ROI MAT files!']);
      end
      clear nroi nmroi;
%
% Loop through T1rho MAT Files
%
      for km = 1:nrho
%
% Load Data
%
         rhonam = rhonams{km};
         load(fullfile(rdir,rhonam),'iszs','nslt','scmx','sns', ...
              'snt','splt','st');
         npix = prod(iszs);  % Number of pixels in an image
         fs = ['S' snt];     % Series number prefaced with a 'S'
%
         idm = contains(roinams,rhonam(1:end-4));     % Get matching file
         roinam = roinams{idm};
         load(fullfile(rdir,roinam),'icmprt','maskf','maskt','rsl');
%
         idm = contains(mroinams,rhonam(1:end-4));    % Get matching file
         mroinam = mroinams{idm};
         load(fullfile(rdir,mroinam),'rsll','rslm');  % Tread slices
%
% Parse Series Text for Leg and Load
%
         if strcmpi(st(1),'L')
           leg = 'L';
           legtxt = 'Left Leg ';
           ileg = 0;    % Coding for leg
         else
           leg = 'R';
           legtxt = 'Right Leg ';
           ileg = 1;
         end
%
         if contains(st,'Load','IgnoreCase',true)
           ld = 'LD';
           ldtxt = 'Loaded';
           ild = 1;     % Coding for load
         else
           ld = 'UL';
           ldtxt = 'Unloaded';
           ild = 0;
         end
%
% Combine Subject, Visit, Leg, Load, and MRI Series for Print Titles
%
         stitle = [subjnam ', ' vnam ', ' legtxt ', ' ldtxt ...
                   ', Series ' snt];
%
% Combine Region Slices and Compartment Identifier
%
         rslb = sort([rsll; rslm]);    % Both compartments
         rsls = {rsll; rslm; rslb};    % Lateral - row 1, medial - row 2, both - row 3
%          nrsls = cellfun('length',rsls);
%          nrsls = cellfun('size',rsls,1);
         nrsls = [size(rsll,1); size(rslm,1); size(rslb,1)];    % Faster than cellfun
%
         [~,idl] = intersect(rsl,rsll);     % Index to lateral compartment
         [~,idm] = intersect(rsl,rslm);     % Index to medial compartment
%
         nidl = size(idl,1);           % Get number of lateral slices
         nidm = size(idm,1);           % Get number of medial slices
%
% Write Out Number of Slices
%
         fprintf(1,'%s, %i, %i\n', stitle, nidl,nidm);
%
      end               % End of km loop - T1rho MAT file loop
%
      fprintf(1,'\n');
%
   end                  % End of kv loop - visits loop
%
   fprintf(1,'\n');
%
end                     % End of ks loop - subjects loop
%
