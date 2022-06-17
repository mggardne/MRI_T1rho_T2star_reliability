%#######################################################################
%
%                     * Combine SHeeTS 2 Program *
%
%          M-File which reads the subject data from different sheets in
%     the MS-Excel spreadsheet file, Registration_data2.xlsx, and
%     combines the subject data into one sheet in the output MS-Excel
%     spreadsheet file, Combined_Registration3.xlsx.
%
%     NOTES:  None.
%
%     22-Apr-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Input and Output MS-Excel Spreadsheet Files
%
ddir = fullfile('Results','NACOB_Final');    % Data directory
xlsf = 'Registration_data2.xlsx';       % Input spreadsheet file name
xlsf = fullfile(ddir,xlsf);
%
xlsf2 = 'Combined_Registration3.xlsx';  % Output spreadsheet file name
xlsf2 = fullfile(ddir,xlsf2);
%
% Get Sheets in Input Spreadsheet File
%
[~,shts] = xlsfinfo(xlsf);
nsht = length(shts);
%
% Loop through Sheets
%
for k = 1:nsht
%
   id = isspace(shts{k});              % Find spaces in sheet name
   id = find(id);
   id = id(1)-1;
   subj = shts{k}(1:id);               % Get subject ID
   v = str2double(shts{k}(end));       % Get visit number
%
% Read Subject Data from Input Spreadsheet File
%
   raw = readcell(xlsf,'Sheet',shts{k});    % Read subject data
%
   nr = size(raw,1);
   raw = [cellstr(repmat(subj,nr,1)) num2cell(repmat(v,nr,1)) raw];
   if k>1
     raw = raw(2:end,:);
   end
%
% Write Subject Data to New Output Spreadsheet File
%
   writecell(raw,xlsf2,'Sheet','Combined','WriteMode','append');
%
end