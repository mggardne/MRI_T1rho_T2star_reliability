%#######################################################################
%
%              * MRI FIT Reliability PTOA Slice Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     MAT files and fits a monoexponential to the MRI data as a function
%     of spin lock or echo times where T1rho or T2* are the time
%     constants of the fits.  Resulting T1rho and T2* values and summary
%     statistics are written to the MS-Excel spreadsheet,
%     mri_fitrp.xlsx, in the "Results\PTOA_ROIs" directory.
%
%     NOTES:  1.  Data MAT files must be in subject directories starting
%             with "MRIR ".
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "prois".  See rd_dicom.m and
%             seg_prois_r.m.
%
%             3.  M-file exp_fun1.m, psl_ana.m and psl_plt.m must be in
%             the current directory or path.
%
%             4.  Currently, set to just do T2* analysis with results
%             in the MS-Excel spreadsheet, mri_fitrp-t2s.xlsx.
%
%     23-Sep-2024 * Mack Gardner-Morse
%
%     16-Apr-2025 * Mack Gardner-Morse * Modified for T2* analysis.
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
resdir = fullfile('Results','PTOA_ROIs');        % Results directory
%
ifirst = true;          % First write to file
% xlsnam = 'mri_fitrp.xlsx';             % Results spreadsheet
xlsnam = 'mri_fitrp_t2s.xlsx';         % Results spreadsheet
xlsnam = fullfile(resdir,xlsnam);      % Include output directory
hdrs1 = {'Subject' 'Visit' 'Result' 'Leg' 'Bone' 'Comprt' 'ROI' ...
         'Layer'};
hdrs2 = {'Pixels' 'T1R/T2S' 'RSS' 'ValidPix' 'Mean' 'Min' 'Max' ...
         'SD' 'COV'};
%
psnam = fullfile(resdir,'mri_fitrp_'); % Start of PS file name
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
% Initialize Results Variables
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit - 1 = Visit 1 and 2 = Visit 2
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Bone - 1 = femur and 2 = tibia
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - ROI - 1 = anterior/trochlea, 2 = central and 3 - posterior
%   Index 7 - Layer - 1 = deep and 2 = superficial
%
t1r_res = zeros(nsubj,2,2,2,2,3,2);
t1r_npx = zeros(nsubj,2,2,2,2,3,2);
t1r_rss = zeros(nsubj,2,2,2,2,3,2);
%
t1r_respx = cell(nsubj,2,2,2,2,3,2);
t1r_rsspx = cell(nsubj,2,2,2,2,3,2);
t1r_nps = cell(nsubj,2,2,2,2,3,2);
%
t2s_res = zeros(nsubj,2,2,2,2,3,2);
t2s_npx = zeros(nsubj,2,2,2,2,3,2);
t2s_rss = zeros(nsubj,2,2,2,2,3,2);
%
t2s_respx = cell(nsubj,2,2,2,2,3,2);
t2s_rsspx = cell(nsubj,2,2,2,2,3,2);
t2s_nps = cell(nsubj,2,2,2,2,3,2);
%
% Loop through Subjects
%
for ks = 1:nsubj
%
% Get Subject Directory (Name) and Number
%
   sdir = sdirs{ks};    % Current subject directory (name)
   subj = eval(sdir(6:7));   % Subject number
   subjtxt = ['Subject ' sdir];
%
   psnams = [psnam sdir(6:7)];         % Add subject to PS file name
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
      vid = kv-1;                      % Visit coding: 0 - Visit 1 and 1 - Visit 2
%
      psnamv = [psnams '_V' vstr];     % Add visit to PS file name
%
% Directory with Data Matlab MAT Files
%
      rdir = fullfile(sdir,vdir);      % Directory with data
%
% Get T1rho MAT Files in Directory
%
   ido = false;         % Skip T1rho
%    ido = true;          % Do T1rho
%
   if ido
%
     d = dir(fullfile(rdir,'T1rho_S*.mat'));
     roinams = {d.name}';
     idr = contains(roinams,'roi','IgnoreCase',true);      % Masks
%
     rhonams = roinams(~idr);
     idx = contains(rhonams,'chk','IgnoreCase',true);      % Check files
     rhonams = rhonams(~idx);          % Image MAT files
     idv = contains(rhonams,'_3');     % Ignore test MAT files
     rhonams = rhonams(~idv);
%
     idr = contains(roinams,'prois.mat','IgnoreCase',true);     % PTOA ROIs
     roinams = roinams(idr);          % PTOA ROI MAT files
     nroi = size(roinams,1);
%
% Get Matching T1rho MAT Files
%
     [~,~,idi] = intersect(regexprep(roinams,'_prois',''),rhonams);
     rhonams = rhonams(idi);
     nrho = size(rhonams,1);
%
     if nrho~=nroi
       error([' *** ERROR in mri_fitrp:  Number of T1rho MAT files', ...
              ' does not match the number of ROI MAT files!']);
     end
     clear nroi;
%
% T1rho Identifier
%
     ires = 0;          % ires = 0 - T1rho, ires = 1 - T2*
     restxt = [subjtxt ' ' vnam ' T1\rho'];
     idt = 1;           % Spin lock/echo time for plots - 1 = 0 ms spin lock time
%
     psnamr = [psnamv '_T1R_'];        % Add result type to PS file name
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
        npix = prod(iszs);   % Number of pixels in an image
        fs = ['S' snt];      % Series number prefaced with a 'S'
%
        idm = contains(roinams,rhonam(1:end-4)); % Get matching file
        roinam = roinams{idm};
        load(fullfile(rdir,roinam),'leg','maskf','maskfr', ...
             'maskt','masktr','nrsl','rsl','rslf','rslt');
%
% Get Leg Text and Coding
%
        if strcmpi(leg,'L')
          legtxt = ' Left Leg';
          ileg = 0;     % Coding for left leg
        else
          legtxt = ' Right Leg';
          ileg = 1;     % Coding for right leg
        end
%
% Add Leg to PS File Name
%
        psnamf = [psnamr leg pstyp];   % Add leg to PS file name
%
% Get Femur ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral-medial, posterior-center, lateral trochlea, and medial trochlea)
%   Third:   Sides of plane (lateral/medial, center/posterior, lateral trochlea/center, and medial trochlea/center)
%   Fourth:  Number of slices
%
        fmasklt = squeeze(maskfr(:,1,1,:)&maskfr(:,3,1,:));     % Lateral trochlea
        fmasklc = squeeze(maskfr(:,1,1,:)&maskfr(:,2,1,:)& ...
                         maskfr(:,3,2,:));       % Lateral central
        fmasklp = squeeze(maskfr(:,1,1,:)&maskfr(:,2,2,:));     % Lateral posterior
%
        fmaskl = {fmasklt; fmasklc; fmasklp};    % Lateral femur
%
        fmaskmt = squeeze(maskfr(:,1,2,:)&maskfr(:,4,1,:));     % Medial trochlea
        fmaskmc = squeeze(maskfr(:,1,2,:)&maskfr(:,2,1,:)& ...
                         maskfr(:,4,2,:));       % Medial central
        fmaskmp = squeeze(maskfr(:,1,2,:)&maskfr(:,2,2,:));     % Medial posterior
%
        fmaskm = {fmaskmt; fmaskmc; fmaskmp};    % Medial femur
%
% Get Tibia ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral central-anterior, lateral posterior-central, medial central-anterior, medial posterior-central)
%   Third:   Sides of plane (anterior/central, or central/posterior)
%   Fourth:  Number of slices
%
        tmaskla = squeeze(masktr(:,1,1,:));      % Lateral anterior
        tmasklc = squeeze(masktr(:,1,2,:)&masktr(:,2,1,:));     % Lateral central
        tmasklp = squeeze(masktr(:,2,2,:));      % Lateral posterior
%
        tmaskl = {tmaskla; tmasklc; tmasklp};    % Lateral tibia
%
        tmaskma = squeeze(masktr(:,3,1,:));      % Medial anterior
        tmaskmc = squeeze(masktr(:,3,2,:)&masktr(:,4,1,:));     % Medial central
        tmaskmp = squeeze(masktr(:,4,2,:));      % Medial posterior
%
        tmaskm = {tmaskma; tmaskmc; tmaskmp};    % Medial tibia
%
% Combine Masks into Cell Arrays
%
        masklay = {maskf; maskt};      % Combine femur and tibia masks
%
        fmask = {fmaskl; fmaskm};      % Combine femur masks
        tmask = {tmaskl; tmaskm};      % Combine tibia masks
        maskroi = {fmask; tmask};      % Combine femur and tibia masks
%
        rslbs = {rslf; rslt};          % Combine femur and tibia slices
%
% Do Slice Analysis
%
        [tc,~,rss,npx,id,tcp,~,rssp,nps] = psl_ana(v,masklay, ...
                     maskroi,rsl,nrsl,rslbs,splt,nslt,fun,init,tr0,opt);
        na = size(tc,1);               % Number of results
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit - 1 = Visit 1 and 2 = Visit 2
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Bone - 1 = femur and 2 = tibia
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - ROI - 1 = anterior/trochlea, 2 = central and 3 - posterior
%   Index 7 - Layer - 1 = deep and 2 = superficial
%
% Note:  Layers for masks and analysis variables are:
%        1 = superficial and 2 = deep.
%
        for ka = 1:na
           t1r_res(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                   id(ka,4)+1) = tc(ka);
           t1r_npx(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                   id(ka,4)+1) = npx(ka);
           t1r_rss(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                   id(ka,4)+1) = rss(ka);
           t1r_respx{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                     id(ka,4)+1} = tcp{ka};
           t1r_rsspx{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                     id(ka,4)+1} = rssp{ka};
           t1r_nps{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                     id(ka,4)+1} = nps{ka};
        end
%
% Plot Results
%
        sid = [restxt legtxt];    % Subject, visit, result (T1rho/T2*), and leg
        psl_plt(v,masklay,maskroi,rsl,nrsl,rslbs,idt,tcp,nps,mxtr, ...
                cmap,sid,psnamf);
%
% Get Statistics on Pixel Results
%
        npxv = zeros(na,1);            % Number of valid results
        tcpm = zeros(na,1);            % Mean
        tcpmn = zeros(na,1);           % Minimum
        tcpmx = zeros(na,1);           % Maximum
        tcpsd = zeros(na,1);           % SD
%
        for ka = 1:na
           idv = tcp{ka}>=trmn&tcp{ka}<=trmx;
           npxv(ka) = sum(idv);        % Number of valid results
           if npxv(ka)>0
             tcpv = tcp{ka}(idv);      % Valid T1rho values
             tcpm(ka) = mean(tcpv);    % Mean
             tcpmn(ka) = min(tcpv);    % Minimum
             tcpmx(ka) = max(tcpv);    % Maximum
             tcpsd(ka) = std(tcpv);    % SD
           end
        end
%
        tcpcov = 100*tcpsd./tcpm;      % Coefficient of variation
        tcpcov(isnan(tcpcov)) = 0;     % Catch any NaNs
%
% Combine Identifiers
%
        ids = [subj vid ires ileg];    % MAT file identifiers
        ids = repmat(ids,na,1);
        ids = [ids id];                % All identifiers
%
% Create and Write Table of Results
%
        t1 = array2table(ids,'VariableNames',hdrs1);
        t2 = table(npx,tc,rss,npxv,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
                   'VariableNames',hdrs2);
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
     end                % End of km loop - T1rho MAT file loop
%
     close all;         % Close all plot windows
%
   end                  % End of ido - Skip T1rho?
%
% Get T2* MAT Files in Directory
%
%    ido = false;         % Skip T2*
   ido = true;          % Do T2*
%
   if ido
%
     d = dir(fullfile(rdir,'T2star_S*.mat'));
     roinams = {d.name}';
     idr = contains(roinams,'roi','IgnoreCase',true);      % Masks
%
     starnams = roinams(~idr);
     idx = contains(starnams,'chk','IgnoreCase',true);     % Check files
     starnams = starnams(~idx);        % Image MAT files
     idv = contains(starnams,'_3');    % Ignore test MAT files
     starnams = starnams(~idv);
%
     idr = contains(roinams,'prois.mat','IgnoreCase',true);     % PTOA ROIs
     roinams = roinams(idr);           % ROI MAT files
     idv = contains(roinams,'_3');     % Ignore test MAT files
     roinams = roinams(~idv);
     nroi = size(roinams,1);
%
% Get Matching T2* MAT Files
%
     [~,~,idi] = intersect(regexprep(roinams,'_prois',''),starnams);
     starnams = starnams(idi);
     nstar = size(starnams,1);
%
     if nstar~=nroi
       error([' *** ERROR in mri_fitrp:  Number of T2* MAT files', ...
              ' does not match the number of ROI MAT files!']);
     end
     clear nroi;
%
% T2* Identifier
%
     ires = 1;          % ires = 0 - T1rho, ires = 1 - T2*
     restxt = [subjtxt ' ' vnam ' T2*'];
     idt = 3;           % Spin lock/echo time for plots - 3 = 5 ms echo time
%
     psnamr = [psnamv '_T2S_'];        % Add result type to PS file name
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
        npix = prod(iszs);   % Number of pixels in an image
        fs = ['S' snt];      % Series number prefaced with a 'S'
%
        idm = contains(roinams,starnam(1:end-4));     % Get matching file
        roinam = roinams{idm};
        load(fullfile(rdir,roinam),'leg','maskf','maskfr', ...
             'maskt','masktr','nrsl','rsl','rslf','rslt');
%
% Parse Series Text for Leg
%
        if strcmpi(leg,'L')
          legtxt = ' Left Leg';
          ileg = 0;       % Coding for leg
        else
          legtxt = ' Right Leg';
          ileg = 1;       % Coding for leg
        end
%
% Add Leg to PS File Name
%
        psnamf = [psnamr leg pstyp];   % Add leg to PS file name
%
% Get Femur ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral-medial, posterior-center, lateral trochlea, and medial trochlea)
%   Third:   Sides of plane (lateral/medial, center/posterior, lateral trochlea/center, and medial trochlea/center)
%   Fourth:  Number of slices
%
        fmasklt = squeeze(maskfr(:,1,1,:)&maskfr(:,3,1,:));     % Lateral trochlea
        fmasklc = squeeze(maskfr(:,1,1,:)&maskfr(:,2,1,:)& ...
                         maskfr(:,3,2,:));       % Lateral central
        fmasklp = squeeze(maskfr(:,1,1,:)&maskfr(:,2,2,:));     % Lateral posterior
%
        fmaskl = {fmasklt; fmasklc; fmasklp};    % Lateral femur
%
        fmaskmt = squeeze(maskfr(:,1,2,:)&maskfr(:,4,1,:));     % Medial trochlea
        fmaskmc = squeeze(maskfr(:,1,2,:)&maskfr(:,2,1,:)& ...
                         maskfr(:,4,2,:));       % Medial central
        fmaskmp = squeeze(maskfr(:,1,2,:)&maskfr(:,2,2,:));     % Medial posterior
%
        fmaskm = {fmaskmt; fmaskmc; fmaskmp};    % Medial femur
%
% Get Tibia ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral central-anterior, lateral posterior-central, medial central-anterior, medial posterior-central)
%   Third:   Sides of plane (anterior/central, or central/posterior)
%   Fourth:  Number of slices
%
        tmaskla = squeeze(masktr(:,1,1,:));      % Lateral anterior
        tmasklc = squeeze(masktr(:,1,2,:)&masktr(:,2,1,:));     % Lateral central
        tmasklp = squeeze(masktr(:,2,2,:));      % Lateral posterior
%
        tmaskl = {tmaskla; tmasklc; tmasklp};    % Lateral tibia
%
        tmaskma = squeeze(masktr(:,3,1,:));      % Medial anterior
        tmaskmc = squeeze(masktr(:,3,2,:)&masktr(:,4,1,:));     % Medial central
        tmaskmp = squeeze(masktr(:,4,2,:));      % Medial posterior
%
        tmaskm = {tmaskma; tmaskmc; tmaskmp};    % Medial tibia
%
% Combine Masks into Cell Arrays
%
        masklay = {maskf; maskt};      % Combine femur and tibia masks
%
        fmask = {fmaskl; fmaskm};      % Combine femur masks
        tmask = {tmaskl; tmaskm};      % Combine tibia masks
        maskroi = {fmask; tmask};      % Combine femur and tibia masks
%
        rslbs = {rslf; rslt};          % Combine femur and tibia slices
%
% Do Slice Analysis
%
        [tc,~,rss,npx,id,tcp,~,rssp,nps] = psl_ana(v,masklay, ...
                       maskroi,rsl,nrsl,rslbs,etns,netn,fun,init,ts0,opt);
        na = size(tc,1);               % Number of results
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit - 1 = Visit 1 and 2 = Visit 2
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Bone - 1 = femur and 2 = tibia
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - ROI - 1 = anterior/trochlea, 2 = central and 3 - posterior
%   Index 7 - Layer - 1 = deep and 2 = superficial
%
% Note:  Layers for masks and compartment analysis variables are:
%        1 = superficial and 2 = deep.
%
        for ka = 1:na
           t2s_res(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                   id(ka,4)+1) = tc(ka);
           t2s_npx(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                   id(ka,4)+1) = npx(ka);
           t2s_rss(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                   id(ka,4)+1) = rss(ka);
           t2s_respx{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                     id(ka,4)+1} = tcp{ka};
           t2s_rsspx{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                     id(ka,4)+1} = rssp{ka};
           t2s_nps{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                     id(ka,4)+1} = nps{ka};
        end
%
% Plot Results
%
        sid = [restxt legtxt];    % Subject, visit, result (T1rho/T2*), and leg
        psl_plt(v,masklay,maskroi,rsl,nrsl,rslbs,idt,tcp,nps,mxts, ...
                cmap,sid,psnamf);
%
% Get Statistics on Pixel Results
%
        npxv = zeros(na,1);            % Number of valid results
        tcpm = zeros(na,1);            % Mean
        tcpmn = zeros(na,1);           % Minimum
        tcpmx = zeros(na,1);           % Maximum
        tcpsd = zeros(na,1);           % SD
%
        for ka = 1:na
           idv = tcp{ka}>=tsmn&tcp{ka}<=tsmx;
           npxv(ka) = sum(idv);        % Number of valid results
           if npxv(ka)>0
             tcpv = tcp{ka}(idv);      % Valid T2* values
             tcpm(ka) = mean(tcpv);    % Mean
             tcpmn(ka) = min(tcpv);    % Minimum
             tcpmx(ka) = max(tcpv);    % Maximum
             tcpsd(ka) = std(tcpv);    % SD
           end
        end
%
        tcpcov = 100*tcpsd./tcpm;      % Coefficient of variation
        tcpcov(isnan(tcpcov)) = 0;     % Catch any NaNs
%
% Combine Identifiers
%
        ids = [subj vid ires ileg];    % MAT file identifiers
        ids = repmat(ids,na,1);
        ids = [ids id];                % All identifiers
%
% Create and Write Table of Results
%
        t1 = array2table(ids,'VariableNames',hdrs1);
        t2 = table(npx,tc,rss,npxv,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
                   'VariableNames',hdrs2);
        t = [t1 t2];
%
%         writetable(t,xlsnam,'WriteMode','append', ...
%                    'WriteVariableNames',false);
        if ifirst
          writetable(t,xlsnam,'WriteMode','replacefile');
          ifirst = false;
        else
          writetable(t,xlsnam,'WriteMode','append', ...
                     'WriteVariableNames',false);
        end
%
     end                % End of km loop - T2* MAT file loop
%
     close all;         % Close all plot windows
%
   end                  % End of ido - Skip T2*?
%
   end                  % End of kv loop - visits loop
%
end                     % End of ks loop - subjects loop
%
% Save to MAT File
%
% save(fullfile(resdir,'mri_fitrp.mat'),'t1r_res','t1r_npx','t1r_rss', ...
%      't1r_respx','t1r_rsspx','t1r_nps','t2s_res','t2s_npx', ...
%      't2s_rss','t2s_respx','t2s_rsspx','t2s_nps');
save(fullfile(resdir,'mri_fitrp_t2s.mat'),'t1r_res','t1r_npx', ...
     't1r_rss','t1r_respx','t1r_rsspx','t1r_nps','t2s_res', ...
     't2s_npx','t2s_rss','t2s_respx','t2s_rsspx','t2s_nps');
%
return