%#######################################################################
%
%                * Combine REGistrations Data Program *
%
%          M-File which reads the registration data from different
%     subject and visit MS-Excel spreadsheet files and combines them as
%     different sheets in the output MS-Excel spreadsheet file,
%     Registration_data?.xlsx.
%
%     NOTES:  1.  The program should be run from the top level
%             directory for the reliability study:
%             \MRI_Reliability_Study\
%
%     14-Jun-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Output MS-Excel Spreadsheet File
%
ddir = fullfile('Results','NACOB_Final');   % Data directory
% xlsf = 'Registration_data.xlsx';      % Input spreadsheet file name
% xlsf = 'Registration_data2.xlsx';      % Input spreadsheet file name
xlsf = 'Registration_data3.xlsx';      % Input spreadsheet file name
xlsf = fullfile(ddir,xlsf);
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
% Loop Through Subject and Visit Spreadsheets
%
for kk = 1:nsubj
%
   for ll = 1:nsdir
%
      xlsdir = fullfile(rdirs{kk},rsdirs{ll});
      dirstr = split(xlsdir,filesep);
      dirstr = [dirstr{end-1} '_' dirstr{end}];
      dirstr = strrep(dirstr,' ','_');
%
      xlsnam = fullfile(xlsdir,[dirstr '.xlsx']);
%
      shtnam = strrep([dirstr(6:end-1) ' ' dirstr(end)],'_',' ');
%
% Read and Write Sheet
%
      raw = readcell(xlsnam);
      raw = raw([1:13 end-15:end],1:13);
      writecell(raw,xlsf,'Sheet',shtnam);
%
   end
end