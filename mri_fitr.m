%#######################################################################
%
%                   * MRI FIT Reliability Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     MAT files and fits a monoexponential to the MRI data as a function
%     of spin lock or echo times where T1rho or T2* are the time
%     constants of the fits.  Resulting T1rho and T2* values and summary
%     statistics are written to the MS-Excel spreadsheet, mri_fitr.xlsx,
%     in the "Results\NACOB" directory.
%
%     NOTES:  1.  Data MAT files must be in subject directories starting
%             with "MRIR" and visit subdirectories "Visit1" or "Visit2".
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "rois".  See rd_dicom.m and seg_rois.m.
%
%             3.  M-file exp_fun1.m, reg_ana.m and reg_plt.m must be in
%             the current directory or path.
%
%             4.  Region is NOT a tibial compartment.  Variable "irsl"
%             indicates the tibial compartment of each region.  The
%             relationship between region and compartment depends on the
%             leg and the direction the slices were digitized.
%
%     27-Jan-2022 * Mack Gardner-Morse
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
% init = -1;              % Use weighted least squares for starting parameters
% init = 0;               % Use linear least squares for starting parameters
init = 1;               % Use fixed starting parameters
tr0 = 65;               % Initial T1rho estimate in ms
trmx = 100;             % Maximum valid T1rho result
trmn = 0;               % Minimum valid T1rho result
ts0 = 35;               % Initial T2* estimate in ms
tsmx = 100;             % Maximum valid T2* result
tsmn = 0;               % Minimum valid T2* result
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results','HonorsThesis');     % Results directory
%
ifirst = true;          % First write to file
xlsnam = 'mri_fitr.xlsx';              % Results spreadsheet
xlsnam = fullfile(resdir,xlsnam);      % Include output directory
hdrs1 = {'Subject' 'Visit' 'Result' 'Leg' 'Load' 'Comprt' 'Bone' ...
         'Layer'};
hdrs2 = {'Pixels', 'T1R/T2S' 'RSS' 'Mean' 'Min' 'Max' 'SD' 'COV'};
%
psnam = fullfile(resdir,'mri_fitr_');  % Start of PS file name
pstyp = '.ps';          % PS file type
%
% Get Subject Directories and Visit Subdirectories
%
sdirs = dir('MRIR*');
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
nsubj = size(sdirs,1);
%
vdirs = {'Visit1'; 'Visit2'};          % Visit directories
nvisit = size(vdirs,1);
%
o3 = ones(1,3);         % Column index for variable "id"
o5 = ones(1,5);         % Column index for variable "ids"
%
% Initialize Results Variables
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - Bone - 1 = femur and 2 = tibia
%   Index 7 - Layer - 1 = deep and 2 = superficial
%
t1r_res = zeros(nsubj,nvisit,2,2,2,2,2);
t1r_npx = zeros(nsubj,nvisit,2,2,2,2,2);
t1r_rss = zeros(nsubj,nvisit,2,2,2,2,2);
%
t1r_respx = cell(nsubj,nvisit,2,2,2,2,2);
t1r_rsspx = cell(nsubj,nvisit,2,2,2,2,2);
t1r_nps = cell(nsubj,nvisit,2,2,2,2,2);
%
t2s_res = zeros(nsubj,nvisit,2,2,2,2,2);
t2s_npx = zeros(nsubj,nvisit,2,2,2,2,2);
t2s_rss = zeros(nsubj,nvisit,2,2,2,2,2);
%
t2s_respx = cell(nsubj,nvisit,2,2,2,2,2);
t2s_rsspx = cell(nsubj,nvisit,2,2,2,2,2);
t2s_nps = cell(nsubj,nvisit,2,2,2,2,2);
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
      nrho = size(rhonams,1);
%
      roinams = roinams(idr);          % ROI MAT files
      nroi = size(roinams,1);
%
      if nrho~=nroi
        error([' *** ERROR in mri_fitr:  Number of T1rho MAT files', ...
               ' does not match the number of ROI MAT files!']);
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
         npix = prod(iszs);  % Number of pixels in an image
         fs = ['S' snt];     % Series number prefaced with a 'S'
%
         idm = contains(roinams,rhonam(1:end-4));     % Get matching file
         roinam = roinams{idm};
         load(fullfile(rdir,roinam),'icmprt','irsl1','irsl2', ...
              'maskfr1','maskfr2','masktr1','masktr2','rsl','rsl1', ...
              'rsl2');
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
         mask1 = {maskfr1; masktr1};   % Combine femur and tibia masks
         mask2 = {maskfr2; masktr2};   % Combine femur and tibia masks
         mask = {mask1; mask2};        % Combine regional masks
%
% Combine Region Slices and Compartment Identifier
%
         rsls = {rsl1; rsl2};          % Region 1 - row 1, region 2 - row 2
         nrsls = [size(rsl1,1); size(rsl1,1)];
         irsl = [irsl1; irsl2];        % Region 1 - row 1, region 2 - row 2
%
% Do Regional Analysis
%
         [tc,~,rss,npx,id,tcp,ampp,rssp,nps] = reg_ana(v,mask,rsls, ...
                                 nrsls,irsl,splt,nslt,fun,init,tr0,opt);
         na = size(tc,1);              % Number of results
%
% Save Results
%
         for ka = 1:na
            t1r_res(ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                    id(ka,3)+1) = tc(ka);
            t1r_npx(ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                    id(ka,3)+1) = npx(ka);
            t1r_rss(ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                    id(ka,3)+1) = rss(ka);
            t1r_respx{ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                      id(ka,3)+1} = tcp{ka};
            t1r_rsspx{ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                      id(ka,3)+1} = rssp{ka};
            t1r_nps{ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                      id(ka,3)+1} = nps{ka};
         end
%
% Plot Results
%
         mxtc = 70;     % Maximum scale on plots
         sid = ['Subject ' subjnam];
         reg_plt(v,mask,rsls,nrsls,idt,tcp,nps,mxtc,cmap,sid,psnamf);
%
% Get Statistics on Pixel Results
%
         npxv = zeros(8,1);            % Number of valid results
         tcpm = zeros(8,1);            % Mean
         tcpmn = zeros(8,1);           % Minimum
         tcpmx = zeros(8,1);           % Maximum
         tcpsd = zeros(8,1);           % SD
%
         for k = 1:na
            idv = tcp{k}>=trmn&tcp{k}<=trmx;
            npxv(k) = sum(idv);        % Number of valid results
            tcpv = tcp{k}(idv);        % Valid T1rho values
            tcpm(k) = mean(tcpv);      % Mean
            tcpmn(k) = min(tcpv);      % Minimum
            tcpmx(k) = max(tcpv);      % Maximum
            tcpsd(k) = std(tcpv);      % SD
         end
%
         tcpcov = 100*tcpsd./tcpm;     % Coefficient of variation
%
% Combine Identifiers
%
         ids = [subj vid ires ileg ild];         % MAT file identifiers
         ids = repmat(ids,na,1);
         ids = [ids id];               % All identifiers
%
% Create and Write Table of Results
%
         t1 = array2table(ids,'VariableNames',hdrs1);
         t2 = table(npxv,tc,rss,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
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
      end               % End of km loop - T1rho MAT file loop
end                     % End of ido - Skip T1rho?
%
% Get T2* MAT Files in Directory
%
      d = dir(fullfile(rdir,'T2star_S*.mat'));
      roinams = {d.name}';
      idr = contains(roinams,'roi','IgnoreCase',true);     % Masks
%
      starnams = roinams(~idr);        % Image MAT files
      nstar = size(starnams,1);
%
      roinams = roinams(idr);          % ROI MAT files
      nroi = size(roinams,1);
%
      if nstar~=nroi
        error([' *** ERROR in mri_fitr:  Number of T2* MAT files', ...
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
         load(fullfile(rdir,roinam),'icmprt','irsl1','irsl2', ...
              'maskfr1','maskfr2','masktr1','masktr2','rsl','rsl1', ...
              'rsl2');
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
         psnamf = [psnamr leg '_' ld pstyp];     % Add leg and load to PS file name
%
% Combine Masks into a Cell Array
%
         mask1 = {maskfr1; masktr1};   % Combine femur and tibia masks
         mask2 = {maskfr2; masktr2};   % Combine femur and tibia masks
         mask = {mask1; mask2};        % Combine regional masks
%
% Combine Region Slices and Compartment Identifier
%
         rsls = {rsl1; rsl2};          % Region 1 - row 1, region 2 - row 2
         nrsls = [size(rsl1,1); size(rsl1,1)];
         irsl = [irsl1; irsl2];        % Region 1 - row 1, region 2 - row 2
%
% Do Regional Analysis
%
         [tc,~,rss,npx,id,tcp,ampp,rssp,nps] = reg_ana(v,mask,rsls, ...
                                 nrsls,irsl,etns,netn,fun,init,ts0,opt);
         na = size(tc,1);              % Number of results
%
% Save Results
%
         for ka = 1:na
            t2s_res(ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                    id(ka,3)+1) = tc(ka);
            t2s_npx(ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                    id(ka,3)+1) = npx(ka);
            t2s_rss(ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                    id(ka,3)+1) = rss(ka);
            t2s_respx{ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                      id(ka,3)+1} = tcp{ka};
            t2s_rsspx{ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                      id(ka,3)+1} = rssp{ka};
            t2s_nps{ks,kv,ileg+1,ild+1,id(ka,1)+1,id(ka,2)+1, ...
                      id(ka,3)+1} = nps{ka};
         end
%
% Plot Results
%
         mxtc = 65;     % Maximum scale on plots
         sid = ['Subject ' subjnam];
         reg_plt(v,mask,rsls,nrsls,idt,tcp,nps,mxtc,cmap,sid,psnamf);
%
% Get Statistics on Pixel Results
%
         npxv = zeros(8,1);            % Number of valid results
         tcpm = zeros(8,1);            % Mean
         tcpmn = zeros(8,1);           % Minimum
         tcpmx = zeros(8,1);           % Maximum
         tcpsd = zeros(8,1);           % SD
%
         for k = 1:na
            idv = tcp{k}>=tsmn&tcp{k}<=tsmx;
            npxv(k) = sum(idv);        % Number of valid results
            tcpv = tcp{k}(idv);        % Valid T2* values
            tcpm(k) = mean(tcpv);      % Mean
            tcpmn(k) = min(tcpv);      % Minimum
            tcpmx(k) = max(tcpv);      % Maximum
            tcpsd(k) = std(tcpv);      % SD
         end
%
         tcpcov = 100*tcpsd./tcpm;     % Coefficient of variation
%
% Combine Identifiers
%
         ids = [subj vid ires ileg ild];         % MAT file identifiers
         ids = repmat(ids,na,1);
         ids = [ids id];               % All identifiers
%
% Create and Write Table of Results
%
         t1 = array2table(ids,'VariableNames',hdrs1);
         t2 = table(npxv,tc,rss,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
                    'VariableNames',hdrs2);
         t = [t1 t2];
%
         writetable(t,xlsnam,'WriteMode','append', ...
                    'WriteVariableNames',false);
%
      end               % End of km loop - T2* MAT file loop
%
   end                  % End of kv loop - visits loop
%
end                     % End of ks loop - subjects loop
%
% Save to MAT File
%
save(fullfile(resdir,'mri_fitr.mat'),'t1r_res','t1r_npx','t1r_rss', ...
     't1r_respx','t1r_rsspx','t1r_nps','t2s_res','t2s_npx', ...
     't2s_rss','t2s_respx','t2s_rsspx','t2s_nps');
%
return