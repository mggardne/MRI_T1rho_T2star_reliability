%
% Get Subject Directories and Visit Subdirectories
%
rdirs = dir('MRIR*');
rdirs = {rdirs([rdirs.isdir]').name}';
nsubj = size(rdirs,1);
%
rsdirs = {'Visit1'; 'Visit2'};
nsdir = 2;
%
tim = zeros(20,1);      % Time to perform registration
%
for kk = 1:nsubj
%
   cd(rdirs{kk});
%
   for ll = 1:nsdir
      mm = nsdir*kk-nsdir+ll;
      cd(rsdirs{ll});
      tstart = tic;
      pwd
      if kk==6&&ll==2
        subj08_AS_V2 = true;           % Use rd_dicomT2s_mask.m to manually run first T2* series for 08-AS on Visit 2
      else
        subj08_AS_V2 = false;
      end
      rd_dicom_m;
      close all;
      tim(mm) = toc(tstart);
      cd ..
   end
%
   cd ..
%
end
%
tim
%
return