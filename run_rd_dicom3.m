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
for kk = 7:nsubj
%
   cd(rdirs{kk});
%
   for ll = 1:nsdir
      mm = nsdir*kk-nsdir+ll;
      cd(rsdirs{ll});
      tstart = tic;
      pwd
      rd_dicom;
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