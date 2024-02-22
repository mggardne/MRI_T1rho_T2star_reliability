%#######################################################################
%
%                       * Series Times Program *
%
%          M-File which reads the DICOM header information for the
%     times of the MRI series in the MRI Reliability Study directories.
%     The elapsed times from the first series to the T1rho and T2*
%     series are written to the MS-Excel spreadsheet, series_time.xlsx,
%     in the "Results" directory.
%
%     NOTES:  1.  Data MAT files must be in subject directories starting
%             with "MRIR" and visit subdirectories "Visit1" or "Visit2".
%
%             2.  Subject AS-08 Visit 2 Left Unloaded T2* has two first
%             echo time series.  The original one of the first echo time
%             has a different field of view and different resolution.
%             Nic resized the image to match the resolution of the
%             other echo times creating a second first series.  The
%             series time for the second one is a duplicate and needs to
%             be manually removed from the final spreadsheet (line 98 
%             duplicates line 94).
%
%     13-Apr-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Output Directory, Output MS-Excel Spreadsheet File and Output Labels
%
resdir = 'Results';     % Results directory
ifirst = true;          % First write to file
% xlsnam = 'series_time.xlsx';           % Results spreadsheet
xlsnam = 'series_times.xlsx';          % Results spreadsheet
xlsnam = fullfile(resdir,xlsnam);      % Include output directory
hdrs = {'Subject' 'Visit' 'Result' 'Leg' 'Load' 'Off_Time'};
% hdrs = {'Subject' 'Visit' 'Result' 'Leg' 'Load' 'E_Time'};
%
% Clock Time of Initial Off Loading
%
ctim0 = { '9:37'    '11:26'            % MRIR-01
         '10:51'     '8:36'            % MRIR-02
          '9:25'     '9:25'            % MRIR-05
         '13:52'    '13:58'            % MRIR-06
         '11:54'    '11:35'            % MRIR-07
         '13:50'    '12:09'            % MRIR-08
         '10:50'    '12:48'            % MRIR-09
          '8:43'     '8:21'            % MRIR-10
         '13:31'    '13:50'            % MRIR-12
         '14:19'    '14:19'};          % MRIR-13
%
% Get Subject Directories and Visit Subdirectories
%
sdirs = dir('MRIR*');
sdirs = {sdirs([sdirs.isdir]').name}';
nsubj = size(sdirs,1);
%
vdirs = {'Visit1'; 'Visit2'};
nvisit = size(vdirs,1);
%
for ks = 1:nsubj
%
   sdir = sdirs{ks};
%
   subjnam = sdir(6:end);              % Subject name as text
   subj = eval(subjnam(1:2));          % Subject number

%
   for kv = 1:nvisit
%
      vdir = vdirs{kv};
      rdir = fullfile(sdir,vdir);      % Directory with data
%
      load(fullfile(rdir,'dicom_lst2.mat'),'afiles','ddirs', ...
           'nimages','stxt');
%
      ndirs = size(ddirs,1);
%
      etim = zeros(ndirs,1);           % Elapsed times
      ctim = cell(ndirs,1);            % Clock times
%
      tim0 = datetime([ctim0{ks,kv} ':00']);     % Initial clock time
%
% Loop through MRI Series Directories and Get Elapsed Times
%
      for k = 1:ndirs
%
         fnam = afiles{k}{1};
%
         info = dicominfo(fullfile(rdir,ddirs{k},fnam));
%
         tim = info.SeriesTime;
         ctim{k} = [tim(1:2) ':' tim(3:4) ':' tim(5:9)];
         tim = datetime(ctim{k});
%          if k==1
%            tim0 = tim;
%          end
%
         etim(k) = minutes(tim-tim0);
%
      end
%
% Use Series Descriptions to Find T1rho Series
%
      idr = find(contains(stxt,'rho','IgnoreCase',true)&nimages>99);
      nr = size(idr,1);      % Number of T1rho series
%
% Use Series Descriptions to Find T2* Series
%
      ids = contains(stxt,'*')|contains(stxt,'T2','IgnoreCase',true);
%
      ids = [ids; false];
      ids = diff(ids);
      ids = find(ids==1)+1;
%
      ns = size(ids,1);      % Number of sets of T2* series
%
% Loop through T1rho Series
%
      atyp = 0;         % T1rho series
%
      for k = 1:nr
%
         idx = idr(k);  % Index to T1rho series
         st = stxt{idx};
%
% Parse Series Text for Leg and Load
%
         if strcmpi(st(1),'L')
           leg = 0;
         else
           leg = 1;
         end
%
         if contains(st,'Load','IgnoreCase',true)
           ld = 1;
         else
           ld = 0;
         end
%
% Write Times to MS-Excel Spreadsheet
%
         id = [subj kv-1 atyp leg ld];
%
         t1 = array2table([id etim(idx)],'VariableNames',hdrs);
         t2 = table({datestr(tim0,15)},ctim(idx), ...
                    'VariableNames',{'Time0' 'Clock_Time'});
         t = [t1 t2];
%
         if ifirst
           writetable(t,xlsnam,'WriteMode','replacefile');
           ifirst = false;
         else
           writetable(t,xlsnam,'WriteMode','append', ...
                      'WriteVariableNames',false);
         end
%
      end               % End of T1rho series
%
% Loop through T2* Series
%
      atyp = 1;         % T2* series
%
      for k = 1:ns
%
         idx = ids(k);  % Index to T2* series
         st = stxt{idx};
%
% Parse Series Text for Leg and Load
%
         if strcmpi(st(1),'L')
           leg = 0;
         else
           leg = 1;
         end
%
         if contains(st,'Load','IgnoreCase',true)
           ld = 1;
         else
           ld = 0;
         end
%
% Write Times to MS-Excel Spreadsheet
%
         id = [subj kv-1 atyp leg ld];
%
         t1 = array2table([id etim(idx)],'VariableNames',hdrs);
         t2 = table({datestr(tim0,15)},ctim(idx), ...
                    'VariableNames',{'Time0' 'Clock_Time'});
         t = [t1 t2];
%
         writetable(t,xlsnam,'WriteMode','append', ...
                    'WriteVariableNames',false);
%
      end               % End of T2* series
%
   end
%
end
%
return