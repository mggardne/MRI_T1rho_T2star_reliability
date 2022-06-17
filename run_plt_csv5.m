%
% Get Subject Directories
%
rdirs = dir('MRIR*');
rdirs = {rdirs([rdirs.isdir]').name}';
nsubj = size(rdirs,1);
%
% for kdir = 1:nsubj
for kdir = 7:nsubj
%
   cd(rdirs{kdir});
   pwd
%
% Run plt_csv5
%
   plt_csv5;
%
   cd ..
%
end
