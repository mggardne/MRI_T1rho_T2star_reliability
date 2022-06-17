%#######################################################################
%
%                     * CHecK DICOM LiST Program *
%
%          M-File which reads dicom_lst2.mat and checks T1rho and T2*
%     series for the correct spin lock and echo times.
%
%     NOTES:  1.  Matlab MAT file dicom_lst2.mat must be in the current
%             directory or path.
%
%     12-Jan-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Load Image File Data from dicom_lst2.mat
%
load dicom_lst2.mat afiles ddirs etn ets idvr isz nimages pspc sn ...
                    splcktc stxt;
%
% Use Series Descriptions to Find T1rho Series
%
idr = contains(stxt,'rho','IgnoreCase',true)&nimages>99;   % T1rho
nr = sum(idr);
%
ddirr = ddirs(idr);     % Subdirectories for T1rho series
nfiler = nimages(idr);  % Numbers of T1rho files
afiler = afiles(idr);   % T1rho files
iszr = isz(idr,:);      % Image sizes in pixels
snr = sn(idr);          % Series numbers
pspcr = pspc(idr,:);    % Pixel sizes
spltr = splcktc(idr);   % T1rho spin lock times
stxtr = stxt(idr);      % T1rho series
%
% Check Spin Lock Times
%
ichk = ~contains(spltr,'0,10,40,80')';
%
for k = find(ichk)
   slt = extractBetween(stxtr(k),'0 ',textBoundary,'Boundaries', ...
                       'inclusive');
   slt = strrep(slt,'ms','');
   slt = strrep(slt,'_',' ');
   slt = strtrim(slt);
   slt = strrep(slt,' ',',');
   spltr(k) = slt;
end
%
% Get Spin Lock Times as a Matrix
%
sltm = cell(1,nr);
%
for k = 1:nr
   sltm{k} = eval(['[' spltr{k} ']'])';
   nslt = size(sltm{k},1);
   if nslt~=4
     warning(' *** rd_dicom:  Incorrect number of spin lock times!');
     sltm{k}
   end
end
%
sltm = cell2mat(sltm)  % Spin lock times in a matrix
%
% Use Series Descriptions to Find T2* Series
%
ids = contains(stxt,'*')|contains(stxt,'T2','IgnoreCase',true); % T2*
ids = [ids; false];
idds = diff(ids);
idss = find(idds==1)+1;
idse = find(idds==-1);
%
ns = length(idss);      % Number of sets of T2* sequences
nse = length(idse);
if ns~=nse
  error([' *** ERROR in rd_dicom:  Number of starting and ending', ...
         ' indices for sets of T2* series are not equal!']);
end
%
ids = [idss idse];      % Index to start of sets (1st column) and end (2nd column)
netn = diff(ids,1,2)+1; % Number of echo times
%
clear idds idse idss nse;
%
% Check Echo Times
%
if ~all(netn==netn(1))
  warning(' *** rd_dicom:  Number of echo times are not the same!');
  netn
end
%
return