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
for k = 1:nsubj
%
   cd(rdirs{k});
%
   for l = 1:nsdir
      cd(rsdirs{l});
      pwd
      chk_dicom_lst;
      cd ..
   end
%
   cd ..
%
end
