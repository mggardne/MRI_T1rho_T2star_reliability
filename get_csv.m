function csvfile = get_csv(csvdir,leg,ld,iro)
%GET_CSV   Reads a directory for sagittal cartilage segmentation CSV
%          files for a particular leg and axial compression load and
%          returns the CSV file name.
%
%          CSVFILE = GET_CSV(CSVDIR,LEG,LD) given the directory name in
%          the string, CSVDIR, either the character 'L' or 'R' for the
%          left or right leg in LEG, and either 'LD' or 'UL' for loaded
%          or unloaded condition in LD, return the CSV file name,
%          CSVFILE, with the matching leg and load.  CSVFILE is an
%          empty cell array if no matching file is found in the
%          directory, CSVDIR.
%
%          CSVFILE = GET_CSV(CSVDIR,LEG,LD,IRO) if IRO is true, the "RO"
%          CSV file is returned; otherwise, the CSV file without "RO"
%          is returned.  The default is to return the "RO" CSV file.
%
%          NOTES:  1.  CSV file names must match the CSV segmentation
%                  naming convention of the MRI reliability study.
%
%                  2.  Only searches for sagittal cartilage
%                  segmentation files.
%
%                  3.  The directory path usually contains the subject
%                  number, visit number, analysis type (T1rho or T2*)
%                  and bone (femur or tibia).
%
%          28-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  iro = true;
end
%
if (nargin<3)
  error(' *** ERROR in GET_CSV:  Three inputs are required!');
end
%
% Get CSV Files
%
csvfile = ['*_' leg '_SAGAR*' ld '*.csv'];
csvfile = dir(fullfile(csvdir,csvfile));
csvfile = {csvfile.name}';
idx = ~contains(csvfile,'dup','IgnoreCase',true);     % Skip duplicate files
csvfile = csvfile(idx);
%
% Check for Files with the Removal of the Overlap (RO Files)
%
idx = contains(csvfile,'RO');
if any(idx)
  if iro
    csvfile = csvfile(idx);
  else
    csvfile = csvfile(~idx);
  end
end
%
nf = size(csvfile,1);
%
if nf>1
  warning(' *** Warning in GET_CSV:  Unique file name not found!');
  fprintf(1,'  CSV files:\n');
  for k = 1:nf
     fprintf(1,'    %s\n',csvfile{k});
  end
  fprintf(1,'\n');
end
%
csvfile = csvfile{1};
%
return
