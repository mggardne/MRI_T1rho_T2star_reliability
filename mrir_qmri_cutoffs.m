%#######################################################################
%
%                 * MRIR Test of qMRI Cutoffs Program *
%
%          M-File which reads the meniscus qMRI results from MAT files,
%     and calculates the means with different cutoffs.  The results
%     are written to the MS-Excel spreadsheets, mri_cutoffs_t1r.xlsx
%     (for T1rho) and mri_cutoffs_t2s.xlsx (for T2*).
%
%     NOTES:  1.  Must start Matlab with the default directory:
%             T1rho\MRI_Reliability_Study\Results\Meniscus
%
%             2. Must have the following qMRI results MAT files:
%             mri_mer_fit_es.mat, and
%             T2S_06Oct2026\mri_merm1_fit_es.mat
%
%     20-Oct-2026 * Mack Gardner-Morse
%

%#######################################################################
%
% Subject Numbers
%
subjn = [ 1
          2
          5
          6
          7
          8
          9
         10
         12
         13 ];
%
% sn = repmat(subjn,1,2,2,2,3,3);
sn = repmat(subjn,72,1);
%
% Visit Numbers
%
v = [zeros(10,1); ones(10,1)];
vn = repmat(v,36,1);
%
% Leg
%
leg = [zeros(20,1); ones(20,1)];
legn = repmat(leg,18,1);
%
% Load
%
ld = [zeros(40,1); ones(40,1)];
ldn = repmat(ld,9,1);
%
% Compartment
%
c = [zeros(80,1); ones(80,1); 2*ones(80,1)];
cn = repmat(c,3,1);
%
% AP
%
apn = [zeros(240,1); ones(240,1); 2*ones(240,1)];
%
% Identifiers
%
id = [sn vn legn ldn cn apn];
idr = [id(:,1:2) zeros(720,1) id(:,3:end)];
%
% Load T1rho qMRI Results MAT File
%
load mri_mer_fit_es.mat;
%
% Loop Through Cutoffs
%
cutoffs = 50:5:100;
nc = size(cutoffs,2);
o = ones(10,2,2,2,3,3);
%
% Initialize Results Matrix
%
t1r = zeros(720,nc);
t1r = [idr t1r];
%
for k = 1:nc
   cutoff = cutoffs(k)*o;
   cutoff = num2cell(cutoff);
   res = cellfun(@(x, cutoff) mean(x(x < cutoff)),t1r_respx,cutoff);
   t1r(:,k+7) = res(:);   
end
%
idn = isnan(t1r);
t1r(idn) = -1;
%
% Load T2* qMRI Results MAT File
%
load T2S_06Oct2026\mri_merm1_fit_es.mat;
%
% Loop Through Cutoffs
%
cutoffs = 40:5:80;
nc = size(cutoffs,2);
%
% Initialize Results Matrix
%
ids = [id(:,1:2) ones(720,1) id(:,3:end)];
t2s = zeros(720,nc);
t2s = [ids t2s];
%
for k = 1:nc
   cutoff = cutoffs(k)*o;
   cutoff = num2cell(cutoff);
   res = cellfun(@(x, cutoff) mean(x(x < cutoff)),t2s_respx,cutoff);
   t2s(:,k+7) = res(:);   
end
%
idn = isnan(t2s);
t2s(idn) = -1;
%
% Write Results to MS-Excel Spreadsheets
%
writematrix(t1r,'mri_cutoffs_t1r.xlsx');
writematrix(t2s,'mri_cutoffs_t2s.xlsx');
%
return