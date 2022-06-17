%#######################################################################
%
%  * CHecK of 3D and 2D Coordinates in Segmentation CSV Files Program *
%
%          M-File which reads the segmentation CSV files to check if
%     the 3D and 2D coordinates in the femur and tibia CSV segmentation
%     files are similar.  The program only checks the CSV files with
%     "_RO" in the filenames (created by rm_overlap.m).  If the 3D and
%     2D coordinates do not match, the CSV filename is printed to the
%     screen.
%
%     NOTES:  1.  M-files decomp.m and rd_roi6.m must be in the current
%             directory or path.
%
%             2.  Future versions of rm_overlap should incorporate this
%             check into the program.
%
%     29-Apr-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Get Subject Directories and Visit Subdirectories
%
sdirs = dir('MRIR*');
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
nsubj = size(sdirs,1);
%
vdirs = {'Visit1'; 'Visit2'};          % Visit directories
nvisit = size(vdirs,1);
%
% Get Analysis and Bone Subdirectories
%
adirs = ['RHO'; 'T2S'];                % Analysis directories
anams = {'T1\rho'; 'T2*'};
%
bdirs = ['Femur'; 'Tibia'];            % Bone directories
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
% Loop through Analysis Directories
%
      for ka = 1:2
%
% Get Analysis Directory
%
         adir = adirs(ka,:);
         anam = anams{ka};
%
% Loop through Bone Directories
%
         for kb = 1:2
%
% Get Bone Directory
%
            bdir = bdirs(kb,:);
%
% Directory with CSV Data Files
%
            rdirk = fullfile(sdir,vdir,adir,bdir);    % Directory with data
%
% Find RO Cartilage Segmentations
%
            dro = dir(rdirk);
            dro = {dro(~[dro.isdir]').name}';    % File names
            idr = contains(dro,'_RO');
            dro = dro(idr);
            nfiles = size(dro,1);
%
% Read RO Cartilage Segmentations
%
            for kr = 1:nfiles
%
               fnam = fullfile(rdirk,dro{kr});
               roi3 = rd_roi6(fnam);
               roi2 = rd_roi6(fnam,true);
%
% Loop through Compartments
%
               nc = size(roi3,1);
%
               for kc = 1:nc
%
                  dat3 = roi3(kc).data';
                  dat2 = roi2(kc).data';
%
% Get 3D to 2D Transformation
%
                  dat21 = dat2{1};
                  npts1 = size(dat21,1);
                  [s,r,t] = decomp(dat3{1},[dat21 zeros(npts1,1)]);
%
% Compare 2D Coordinates
%
                  dat2 = cell2mat(dat2);
                  dat2t = cell2mat(dat3);
                  npts = size(dat2t,1);
                  dat2t = s*dat2t*r+repmat(t,npts,1);
                  dat2t = dat2t(:,1:2);
                  dd = dat2t-dat2;
                  dd = max(abs(dd(:)));
                  if dd>0.1
                    fprintf(1,' *** 3D/2D Mismatch in file: %s\n',fnam);
                  end                  
%
               end      % End of compartment loop - kc
%
            end         % End of RO CSV file loop - kr
%
         end            % End of bone loop - kb
%
      end               % End of analysis loop - ka
%
   end                  % End of visit loop - kv
%
end                     % End of subject loop - ks
%
return