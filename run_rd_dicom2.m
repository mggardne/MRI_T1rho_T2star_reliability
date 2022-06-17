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
tim = zeros(4,1);       % Time to perform registration
%
for kk = [3 4]
%
   cd(rdirs{k});
%
   if kk==1
     ll = 1;
   else
     ll = 2;
   end
%
   cd(rsdirs{ll});
   tstart = tic;
   pwd
   rd_dicom;
   close all;
   tim(kk) = toc(tstart);
   cd ..
%
   cd ..
%
end
%
tim
%
return