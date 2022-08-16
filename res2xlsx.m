%#######################################################################
%
%           * RESults MAT File to XLSX Spreadsheet Program *
%
%          M-File which reads the MRI reliability results MAT file and
%     outputs the results to the output MS-Excel spreadsheet file,
%     mri_fitr3mat.xlsx.
%
%     NOTES:  1.  The program should be run from the results directory
%             for the reliability study:
%             \MRI_Reliability_Study\Results\NACOB_Final
%
%     19-Jul-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Get Results MAT File Data
%
load mri_fitr3.mat;
%
% MS-Excel Spreadsheet and Column Headers
%
xlsnam = 'mri_fitr3mat.xlsx';          % Spreadsheet file name
hdrs1 = {'Subject' 'Visit' 'Result' 'Leg' 'Load' 'Comprt' 'Bone' ...
'Layer'};
hdrs2 = {'Pixels' 'T1R/T2S' 'RSS' 'ValidPix' 'Mean' 'Min' 'Max' ...
'SD' 'COV'};
%
% T1rho Data from Fit of Entire ROIs
%
dat = [t1r_npx(:) t1r_res(:) t1r_rss(:)];   % Fit over ROI
%
% Get T1rho Statistics on Pixel Results
%
id0 = zeros(640,8);
isz = size(t1r_respx);
trmn = 0;
trmx = 100;
npxv = zeros(640,1);    % Number of valid results
tcpm = zeros(640,1);    % Mean
tcpmn = zeros(640,1);   % Minimum
tcpmx = zeros(640,1);   % Maximum
tcpsd = zeros(640,1);   % SD
%
for k = 1:640
   [i1,i2,i3,i4,i5,i6,i7] = ind2sub(isz,k);
   id0(k,[1:2 4:8]) = [i1,i2,i3,i4,i5,i6,i7];
   res = t1r_respx{i1,i2,i3,i4,i5,i6,i7};
   idv = res>=trmn&res<=trmx;
   npxv(k) = sum(idv);                 % Number of valid results
   tcpv = res(idv);                    % Valid T1rho values
   tcpm(k) = mean(tcpv);               % Mean
   tcpmn(k) = min(tcpv);               % Minimum
   tcpmx(k) = max(tcpv);               % Maximum
   tcpsd(k) = std(tcpv);               % SD
end
%
id0(:,2) = id0(:,2)-1;
id0(:,4:8) = id0(:,4:8)-1;
%
tcpcov = 100*tcpsd./tcpm;              % Coefficient of variation
%
% Create and Write Table of Results to Spreadsheet
%
t1 = array2table(id0,'VariableNames',hdrs1);
t2 = table(dat(:,1),dat(:,2),dat(:,3),npxv,tcpm,tcpmn,tcpmx,tcpsd, ...
           tcpcov,'VariableNames',hdrs2);
t = [t1 t2];
writetable(t,xlsnam,'WriteMode','replacefile');
%
% T2* Data from Fit of Entire ROIs
%
dat = [t2s_npx(:) t2s_res(:) t2s_rss(:)];   % Fit over ROI
%
% Get T2* Statistics on Pixel Results
%
id0 = zeros(640,8);
id0(:,3) = 1;           % T2* analysis
isz = size(t2s_respx);
trmn = 0;
trmx = 100;
npxv = zeros(640,1);    % Number of valid results
tcpm = zeros(640,1);    % Mean
tcpmn = zeros(640,1);   % Minimum
tcpmx = zeros(640,1);   % Maximum
tcpsd = zeros(640,1);   % SD
%
for k = 1:640
   [i1,i2,i3,i4,i5,i6,i7] = ind2sub(isz,k);
   id0(k,[1:2 4:8]) = [i1,i2,i3,i4,i5,i6,i7];
   res = t2s_respx{i1,i2,i3,i4,i5,i6,i7};
   idv = res>=trmn&res<=trmx;
   npxv(k) = sum(idv);                 % Number of valid results
   tcpv = res(idv);                    % Valid T2* values
   tcpm(k) = mean(tcpv);               % Mean
   tcpmn(k) = min(tcpv);               % Minimum
   tcpmx(k) = max(tcpv);               % Maximum
   tcpsd(k) = std(tcpv);               % SD
end
%
id0(:,2) = id0(:,2)-1;
id0(:,4:8) = id0(:,4:8)-1;
%
tcpcov = 100*tcpsd./tcpm;              % Coefficient of variation
%
% Write Table of Results to Spreadsheet
%
t1 = array2table(id0,'VariableNames',hdrs1);
t2 = table(dat(:,1),dat(:,2),dat(:,3),npxv,tcpm,tcpmn,tcpmx,tcpsd, ...
           tcpcov,'VariableNames',hdrs2);
t = [t1 t2];
writetable(t,xlsnam,'WriteMode','append','WriteVariableNames',false);
%
return