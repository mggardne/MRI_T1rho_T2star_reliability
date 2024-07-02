%#######################################################################
%
%             * MRI FIT Reliability 4 Unloaded 2 Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     MAT files and fits a monoexponential to the unloaded MRI data as
%     a function of spin lock or echo times where T1rho or T2* are the
%     time constants of the fits.  Resulting T1rho and T2* values and
%     summary statistics are written to the MS-Excel spreadsheets,
%     mri_fitr4u21.xlsx (loaded ROIs) and mri_fitr4u22.xlsx (unloaded
%     femur ROIs), in the "Results\Unloaded_Femur4u2" directory.
%
%          The unloaded (noncontact) cartilage region of interests are
%     regions on the posterior femoral condyles.  See seg_2rois2.m.
%
%          This analysis combines the deep and superficial cartilage
%     layers. 
%
%     NOTES:  1.  Data MAT files must be in subject directories starting
%             with "MRIR " and visit subdirectories "Visit1" or
%             "Visit2".
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "2rois2".  See rd_dicom.m and
%             seg_2rois2.m.
%
%             3.  M-file exp_fun1.m, cmprt_ana4.m, cmprt_ana4u.m, and
%             cmprt_plt4u.m must be in the current directory or path.
%
%             4.  Update of mri_fitr4u.m to read the new segmentation
%             MAT files with fewer tibial eminence slices. 
%
%     07-Jun-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Setup Color Map
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
%
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter', ...
               2e+3,'Algorithm','levenberg-marquardt','Jacobian', ...
               'on','UseParallel',true);
%
fun = @exp_fun1;        % Exponential function
%
% Initialize Parameters
%
% init = -1;              % Use weighted least squares for starting parameters
% init = 0;               % Use linear least squares for starting parameters
init = 1;               % Use fixed starting parameters
tr0 = 65;               % Initial T1rho estimate in ms
% tr0 = 80;               % Initial T1rho estimate in ms
trmx = 100;             % Maximum valid T1rho result
trmn = 0;               % Minimum valid T1rho result
ts0 = 35;               % Initial T2* estimate in ms
tsmx = 100;             % Maximum valid T2* result
tsmn = 0;               % Minimum valid T2* result
%
mxtr = 80;              % Maximum scale on T1rho plots
mxts = 75;              % Maximum scale on T2* plots
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results','Unloaded_Femur4u2');     % Results directory
%
ifirst1 = true;         % First write to file
xlsnam1 = 'mri_fitr4u21.xlsx';         % Loaded results spreadsheet
xlsnam1 = fullfile(resdir,xlsnam1);    % Include output directory
%
ifirst2 = true;         % First write to file
xlsnam2 = 'mri_fitr4u22.xlsx';         % Unloaded results spreadsheet
xlsnam2 = fullfile(resdir,xlsnam2);    % Include output directory
%
hdrs1_1 = {'Subject' 'Visit' 'Result' 'Leg' 'Load' 'Comprt' 'Bone'};
hdrs1_2 = {'Pixels' 'T1R/T2S' 'RSS' 'ValidPix' 'Mean' 'Min' 'Max' ...
           'SD' 'COV'};
%
hdrs2_1 = {'Subject' 'Visit' 'Result' 'Leg' 'Load' 'Comprt'};
%
psnam = fullfile(resdir,'mri_fitr4u2_');    % Start of PS file name
pstyp = '.ps';          % PS file type
%
% Get Subject Directories and Visit Subdirectories
%
sdirs = dir('MRIR *');
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
nsubj = size(sdirs,1);
%
vdirs = {'Visit1'; 'Visit2'};          % Visit directories
nvisit = size(vdirs,1);
%
% Initialize Loaded Results Variables
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - Bone - 1 = femur and 2 = tibia
%
t1r_res = zeros(nsubj,nvisit,2,2,2,2);
t1r_npx = zeros(nsubj,nvisit,2,2,2,2);
t1r_rss = zeros(nsubj,nvisit,2,2,2,2);
%
t1r_respx = cell(nsubj,nvisit,2,2,2,2);
t1r_rsspx = cell(nsubj,nvisit,2,2,2,2);
t1r_nps = cell(nsubj,nvisit,2,2,2,2);
%
t2s_res = zeros(nsubj,nvisit,2,2,2,2);
t2s_npx = zeros(nsubj,nvisit,2,2,2,2);
t2s_rss = zeros(nsubj,nvisit,2,2,2,2);
%
t2s_respx = cell(nsubj,nvisit,2,2,2,2);
t2s_rsspx = cell(nsubj,nvisit,2,2,2,2);
t2s_nps = cell(nsubj,nvisit,2,2,2,2);
%
% Initialize Unloaded Results Variables
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%
t1ru_res = zeros(nsubj,nvisit,2,2,2);
t1ru_npx = zeros(nsubj,nvisit,2,2,2);
t1ru_rss = zeros(nsubj,nvisit,2,2,2);
%
t1ru_respx = cell(nsubj,nvisit,2,2,2);
t1ru_rsspx = cell(nsubj,nvisit,2,2,2);
t1ru_nps = cell(nsubj,nvisit,2,2,2);
%
t2su_res = zeros(nsubj,nvisit,2,2,2);
t2su_npx = zeros(nsubj,nvisit,2,2,2);
t2su_rss = zeros(nsubj,nvisit,2,2,2);
%
t2su_respx = cell(nsubj,nvisit,2,2,2);
t2su_rsspx = cell(nsubj,nvisit,2,2,2);
t2su_nps = cell(nsubj,nvisit,2,2,2);
%
% Loop through Subjects
%
for ks = 1:nsubj
%
% Get Subject Directory, Name and Number
%
   sdir = sdirs{ks};                   % Current subject directory
   subjnam = sdir(6:end);              % Subject name as text
   subj = eval(subjnam(1:2));          % Subject number
%
   psnams = [psnam subjnam];           % Add subject to PS file name
%
% Loop through Visits
%
   for kv = 1:nvisit
%
% Get Visit Subdirectory, Name and Number
%
      vdir = vdirs{kv};                % Current visit directory
      vstr = int2str(kv);              % Visit number as text
      vnam = ['Visit ' vstr];          % Visit name as text
      vid = kv-1;                      % Visit number
%
      psnamv = [psnams '_V' vstr];     % Add visit to PS file name
%
% Directory with Data Matlab MAT Files
%
      rdir = fullfile(sdir,vdir);      % Directory with data
%
% Get T1rho MAT Files in Directory
%
% ido = false;            % Skip T1rho
ido = true;             % Do T1rho
%
if ido
%
      d = dir(fullfile(rdir,'T1rho_S*.mat'));
      roinams = {d.name}';
      idr = contains(roinams,'roi','IgnoreCase',true);     % Masks
%
      rhonams = roinams(~idr);         % Image MAT files
      idv = ~contains(rhonams,'_3');   % Ignore test MAT files
      rhonams = rhonams(idv);
      nrho = size(rhonams,1);
%
      idr = contains(roinams,'2rois2.mat','IgnoreCase',true);   % Small ROIs
      roinams = roinams(idr);          % ROI MAT files
      nroi = size(roinams,1);
%
      if nrho~=nroi
        error([' *** ERROR in mri_fitr4u2:  Number of T1rho MAT', ...
               ' files does not match the number of ROI MAT files!']);
      end
      clear nroi;
%
% T1rho Identifier
%
      ires = 0;         % ires = 0 - T1rho, ires = 1 - T2*
      idt = 1;          % Spin lock/echo time for plots - 1 = 0 ms spin lock time
%
      psnamr = [psnamv '_T1R_'];       % Add result type to PS file name
%
% Loop through T1rho MAT Files
%
      for km = 1:nrho
%
% Load Data
%
         rhonam = rhonams{km};
         load(fullfile(rdir,rhonam),'iszs','nslt','scmx','sns', ...
              'snt','splt','st','v');
         fs = ['S' snt];     % Series number prefaced with a 'S'
%
         idm = contains(roinams,rhonam(1:end-4));     % Get matching file
         roinam = roinams{idm};
         load(fullfile(rdir,roinam),'maskfrl1','maskfrm1', ...
              'maskfrl2','maskfrm2','masktrl1','masktrm1','rsll1', ...
              'rslm1','rsll2','rslm2');
%
% Parse Series Text for Leg and Load
%
         if strcmpi(st(1),'L')
           leg = 'L';
           ileg = 0;    % Coding for leg
         else
           leg = 'R';
           ileg = 1;
         end
%
         if contains(st,'Load','IgnoreCase',true)
           ld = 'LD';
           ild = 1;     % Coding for load
         else
           ld = 'UL';
           ild = 0;
         end
%
% Add Leg and Load to PS File Name
%
         psnamf = [psnamr leg '_' ld pstyp];     % Add leg and load to PS file name
%
% Combine Masks into a Cell Array
%
         maskl1 = {maskfrl1; masktrl1};     % Combine femur and tibia masks
         maskm1 = {maskfrm1; masktrm1};     % Combine femur and tibia masks
         mask1 = {maskl1; maskm1};          % Combine compartment masks
%
         mask2 = {maskfrl2; maskfrm2};      % Combine compartment masks
%
% Combine Region Slices and Compartment Identifier
%
         rsls1 = {rsll1; rslm1};       % Lateral - row 1, medial - row 2
         nrsls1 = [size(rsll1,1); size(rslm1,1)];
%
         rsls2 = {rsll2; rslm2};       % Lateral - row 1, medial - row 2
         nrsls2 = [size(rsll2,1); size(rslm2,1)];
%
% Do Loaded Compartmental Analysis
%
         [tc1,~,rss1,npx1,id1,tcp1,ampp1,rssp1,nps1] = cmprt_ana4(v, ...
                         mask1,rsls1,nrsls1,splt,nslt,fun,init,tr0,opt);
         na1 = size(tc1,1);            % Number of results
%
% Save Loaded Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - Bone - 1 = femur and 2 = tibia
%
         for ka = 1:na1
            t1r_res(ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1) = ...
                    tc1(ka);
            t1r_npx(ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1) = ...
                    npx1(ka);
            t1r_rss(ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1) = ...
                    rss1(ka);
            t1r_respx{ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1} = ...
                      tcp1{ka};
            t1r_rsspx{ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1} = ...
                      rssp1{ka};
            t1r_nps{ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1} = ...
                    nps1{ka};
         end
%
% Get Statistics on Pixel Results
%
         npxv = zeros(na1,1);          % Number of valid results
         tcpm = zeros(na1,1);          % Mean
         tcpmn = zeros(na1,1);         % Minimum
         tcpmx = zeros(na1,1);         % Maximum
         tcpsd = zeros(na1,1);         % SD
%
         for ka = 1:na1
            idv = tcp1{ka}>=trmn&tcp1{ka}<=trmx;
            npxv(ka) = sum(idv);       % Number of valid results
            tcpv = tcp1{ka}(idv);      % Valid T1rho values
            tcpm(ka) = mean(tcpv);     % Mean
            tcpmn(ka) = min(tcpv);     % Minimum
            tcpmx(ka) = max(tcpv);     % Maximum
            tcpsd(ka) = std(tcpv);     % SD
         end
%
         tcpcov = 100*tcpsd./tcpm;     % Coefficient of variation
%
% Combine Identifiers
%
         ids = [subj vid ires ileg ild];         % MAT file identifiers
         ids = repmat(ids,na1,1);
         ids = [ids id1];               % All identifiers
%
% Create and Write Table of Results
%
         t1 = array2table(ids,'VariableNames',hdrs1_1);
         t2 = table(npx1,tc1,rss1,npxv,tcpm,tcpmn,tcpmx,tcpsd, ...
                    tcpcov,'VariableNames',hdrs1_2);
         t = [t1 t2];
%
         if ifirst1
           writetable(t,xlsnam1,'WriteMode','replacefile');
           ifirst1 = false;
         else
           writetable(t,xlsnam1,'WriteMode','append', ...
                      'WriteVariableNames',false);
         end
%
% Do Unloaded Compartmental Analysis
%
         [tc2,~,rss2,npx2,id2,tcp2,ampp2,rssp2,nps2] = cmprt_ana4u( ...
                       v,mask2,rsls2,nrsls2,splt,nslt,fun,init,tr0,opt);
         na2 = size(tc2,1);            % Number of results
%
% Save Unloaded Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%
         for ka = 1:na2
            t1ru_res(ks,kv,ileg+1,ild+1,id2(ka)+1) = tc2(ka);
            t1ru_npx(ks,kv,ileg+1,ild+1,id2(ka)+1) = npx2(ka);
            t1ru_rss(ks,kv,ileg+1,ild+1,id2(ka)+1) = rss2(ka);
            t1ru_respx{ks,kv,ileg+1,ild+1,id2(ka)+1} = tcp2{ka};
            t1ru_rsspx{ks,kv,ileg+1,ild+1,id2(ka)+1} = rssp2{ka};
            t1ru_nps{ks,kv,ileg+1,ild+1,id2(ka)+1} = nps2{ka};
         end
%
% Get Statistics on Pixel Results
%
         npxv = zeros(na2,1);          % Number of valid results
         tcpm = zeros(na2,1);          % Mean
         tcpmn = zeros(na2,1);         % Minimum
         tcpmx = zeros(na2,1);         % Maximum
         tcpsd = zeros(na2,1);         % SD
%
         for ka = 1:na2
            idv = tcp2{ka}>=trmn&tcp2{ka}<=trmx;
            npxv(ka) = sum(idv);       % Number of valid results
            tcpv = tcp2{ka}(idv);      % Valid T1rho values
            tcpm(ka) = mean(tcpv);     % Mean
            tcpmn(ka) = min(tcpv);     % Minimum
            tcpmx(ka) = max(tcpv);     % Maximum
            tcpsd(ka) = std(tcpv);     % SD
         end
%
         tcpcov = 100*tcpsd./tcpm;     % Coefficient of variation
%
% Combine Identifiers
%
         ids = [subj vid ires ileg ild];         % MAT file identifiers
         ids = repmat(ids,na2,1);
         ids = [ids id2];               % All identifiers
%
% Create and Write Table of Results
%
         t1 = array2table(ids,'VariableNames',hdrs2_1);
         t2 = table(npx2,tc2,rss2,npxv,tcpm,tcpmn,tcpmx,tcpsd, ...
                    tcpcov,'VariableNames',hdrs1_2);
         t = [t1 t2];
%
         if ifirst2
           writetable(t,xlsnam2,'WriteMode','replacefile');
           ifirst2 = false;
         else
           writetable(t,xlsnam2,'WriteMode','append', ...
                      'WriteVariableNames',false);
         end
%
% Plot Results
%
         sid = ['Subject ' subjnam];
         cmprt_plt4u(v,mask1,rsls1,mask2,rsls2,idt,tcp1,nps1,tcp2, ...
                     nps2,mxtr,cmap,sid,psnamf);
%
      end               % End of km loop - T1rho MAT file loop
%
      close all;        % Close all plot windows
%
end                     % End of ido - Skip T1rho?
%
% Get T2* MAT Files in Directory
%
      d = dir(fullfile(rdir,'T2star_S*.mat'));
      roinams = {d.name}';
      idr = contains(roinams,'roi','IgnoreCase',true);     % Masks
%
      starnams = roinams(~idr);        % Image MAT files
      idv = ~contains(starnams,'_3');  % Ignore test MAT files
      starnams = starnams(idv);
      nstar = size(starnams,1);
%
      idr = contains(roinams,'2rois2.mat','IgnoreCase',true);   % Small ROIs
      roinams = roinams(idr);          % ROI MAT files
      nroi = size(roinams,1);
%
      if nstar~=nroi
        error([' *** ERROR in mri_fitr4u2:  Number of T2* MAT files', ...
               ' does not match the number of ROI MAT files!']);
      end
      clear nroi;
%
% T2* Identifier
%
      ires = 1;         % ires = 0 - T1rho, ires = 1 - T2*
      idt = 3;          % Spin lock/echo time for plots - 3 = 5 ms echo time
%
      psnamr = [psnamv '_T2S_'];       % Add result type to PS file name
%
% Loop through T2* MAT Files
%
      for km = 1:nstar
%
% Load Data
%
         starnam = starnams{km};
         load(fullfile(rdir,starnam),'etns','iszs','netn','scmx', ...
              'sns','snt','st','v');
         fs = ['S' snt];     % Series number prefaced with a 'S'
%
         idm = contains(roinams,starnam(1:end-4));    % Get matching file
         roinam = roinams{idm};
         load(fullfile(rdir,roinam),'maskfrl1','maskfrm1', ...
              'maskfrl2','maskfrm2','masktrl1','masktrm1','rsll1', ...
              'rslm1','rsll2','rslm2');
%
% Parse Series Text for Leg and Load
%
         if strcmpi(st(1),'L')
           leg = 'L';
           ileg = 0;    % Coding for leg
         else
           leg = 'R';
           ileg = 1;
         end
%
         if contains(st,'Load','IgnoreCase',true)
           ld = 'LD';
           ild = 1;     % Coding for load
         else
           ld = 'UL';
           ild = 0;
         end
%
% Add Leg and Load to PS File Name
%
         psnamf = [psnamr leg '_' ld pstyp];     % Add leg and load to PS file name
%
% Combine Masks into a Cell Array
%
         maskl1 = {maskfrl1; masktrl1};     % Combine femur and tibia masks
         maskm1 = {maskfrm1; masktrm1};     % Combine femur and tibia masks
         mask1 = {maskl1; maskm1};          % Combine compartment masks
%
         mask2 = {maskfrl2; maskfrm2};      % Combine compartment masks
%
% Combine Region Slices and Compartment Identifier
%
         rsls1 = {rsll1; rslm1};       % Lateral - row 1, medial - row 2
         nrsls1 = [size(rsll1,1); size(rslm1,1)];
%
         rsls2 = {rsll2; rslm2};       % Lateral - row 1, medial - row 2
         nrsls2 = [size(rsll2,1); size(rslm2,1)];
%
% Do Compartmental Analysis
%
         [tc1,~,rss1,npx1,id1,tcp1,ampp1,rssp1,nps1] = cmprt_ana4(v, ...
                         mask1,rsls1,nrsls1,etns,netn,fun,init,ts0,opt);
         na1 = size(tc1,1);              % Number of results
%
% Save Loaded Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - Bone - 1 = femur and 2 = tibia
%
         for ka = 1:na1
            t2s_res(ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1) = ...
                    tc1(ka);
            t2s_npx(ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1) = ...
                    npx1(ka);
            t2s_rss(ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1) = ...
                    rss1(ka);
            t2s_respx{ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1} = ...
                      tcp1{ka};
            t2s_rsspx{ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1} = ...
                      rssp1{ka};
            t2s_nps{ks,kv,ileg+1,ild+1,id1(ka,1)+1,id1(ka,2)+1} = ...
                    nps1{ka};
         end
%
% Get Statistics on Pixel Results
%
         npxv = zeros(na1,1);          % Number of valid results
         tcpm = zeros(na1,1);          % Mean
         tcpmn = zeros(na1,1);         % Minimum
         tcpmx = zeros(na1,1);         % Maximum
         tcpsd = zeros(na1,1);         % SD
%
         for ka = 1:na1
            idv = tcp1{ka}>=tsmn&tcp1{ka}<=tsmx;
            npxv(ka) = sum(idv);       % Number of valid results
            tcpv = tcp1{ka}(idv);      % Valid T2* values
            tcpm(ka) = mean(tcpv);     % Mean
            tcpmn(ka) = min(tcpv);     % Minimum
            tcpmx(ka) = max(tcpv);     % Maximum
            tcpsd(ka) = std(tcpv);     % SD
         end
%
         tcpcov = 100*tcpsd./tcpm;     % Coefficient of variation
%
% Combine Identifiers
%
         ids = [subj vid ires ileg ild];         % MAT file identifiers
         ids = repmat(ids,na1,1);
         ids = [ids id1];              % All identifiers
%
% Create and Write Table of Results
%
         t1 = array2table(ids,'VariableNames',hdrs1_1);
         t2 = table(npx1,tc1,rss1,npxv,tcpm,tcpmn,tcpmx,tcpsd, ...
                    tcpcov,'VariableNames',hdrs1_2);
         t = [t1 t2];
%
         writetable(t,xlsnam1,'WriteMode','append', ...
                    'WriteVariableNames',false);
%
% Do Unloaded Compartmental Analysis
%
         [tc2,~,rss2,npx2,id2,tcp2,ampp2,rssp2,nps2] = cmprt_ana4u( ...
                       v,mask2,rsls2,nrsls2,etns,netn,fun,init,tr0,opt);
         na2 = size(tc2,1);            % Number of results
%
% Save Unloaded Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%
         for ka = 1:na2
            t2su_res(ks,kv,ileg+1,ild+1,id2(ka)+1) = tc2(ka);
            t2su_npx(ks,kv,ileg+1,ild+1,id2(ka)+1) = npx2(ka);
            t2su_rss(ks,kv,ileg+1,ild+1,id2(ka)+1) = rss2(ka);
            t2su_respx{ks,kv,ileg+1,ild+1,id2(ka)+1} = tcp2{ka};
            t2su_rsspx{ks,kv,ileg+1,ild+1,id2(ka)+1} = rssp2{ka};
            t2su_nps{ks,kv,ileg+1,ild+1,id2(ka)+1} = nps2{ka};
         end
%
% Get Statistics on Pixel Results
%
         npxv = zeros(na2,1);          % Number of valid results
         tcpm = zeros(na2,1);          % Mean
         tcpmn = zeros(na2,1);         % Minimum
         tcpmx = zeros(na2,1);         % Maximum
         tcpsd = zeros(na2,1);         % SD
%
         for ka = 1:na2
            idv = tcp2{ka}>=trmn&tcp2{ka}<=trmx;
            npxv(ka) = sum(idv);       % Number of valid results
            tcpv = tcp2{ka}(idv);      % Valid T1rho values
            tcpm(ka) = mean(tcpv);     % Mean
            tcpmn(ka) = min(tcpv);     % Minimum
            tcpmx(ka) = max(tcpv);     % Maximum
            tcpsd(ka) = std(tcpv);     % SD
         end
%
         tcpcov = 100*tcpsd./tcpm;     % Coefficient of variation
%
% Combine Identifiers
%
         ids = [subj vid ires ileg ild];         % MAT file identifiers
         ids = repmat(ids,na2,1);
         ids = [ids id2];               % All identifiers
%
% Create and Write Table of Results
%
         t1 = array2table(ids,'VariableNames',hdrs2_1);
         t2 = table(npx2,tc2,rss2,npxv,tcpm,tcpmn,tcpmx,tcpsd, ...
                    tcpcov,'VariableNames',hdrs1_2);
         t = [t1 t2];
%
         if ifirst2
           writetable(t,xlsnam2,'WriteMode','replacefile');
           ifirst2 = false;
         else
           writetable(t,xlsnam2,'WriteMode','append', ...
                      'WriteVariableNames',false);
         end
%
% Plot Results
%
         sid = ['Subject ' subjnam];
         cmprt_plt4u(v,mask1,rsls1,mask2,rsls2,idt,tcp1,nps1,tcp2, ...
                     nps2,mxtr,cmap,sid,psnamf);
%
      end               % End of km loop - T2* MAT file loop
%
      close all;        % Close all plot windows
%
   end                  % End of kv loop - visits loop
%
end                     % End of ks loop - subjects loop
%
% Save to MAT File
%
save(fullfile(resdir,'mri_fitr4u2.mat'),'t1r_res','t1r_npx', ...
     't1r_rss','t1r_respx','t1r_rsspx','t1r_nps','t2s_res', ...
     't2s_npx','t2s_rss','t2s_respx','t2s_rsspx','t2s_nps', ...
     't1ru_res','t1ru_npx','t1ru_rss','t1ru_respx','t1ru_rsspx', ...
     't1ru_nps','t2su_res','t2su_npx','t2su_rss','t2su_respx', ...
     't2su_rsspx','t2su_nps');
%
return