%#######################################################################
%
%                  * MRI FIT Reliability All Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     MAT files and fits a monoexponential to the MRI data as a function
%     of spin lock or echo times where T1rho or T2* are the time
%     constants of the fits.  Resulting T1rho and T2* values and summary
%     statistics are written to the MS-Excel spreadsheet,
%     mri_fitra.xlsx, in the "Results\All" directory.
%
%     NOTES:  1.  Data MAT files must be in subject directories starting
%             with "MRIR" and visit subdirectories "Visit1" or "Visit2".
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "2rois".  See rd_dicom.m and
%             seg_2rois.m.
%
%             3.  M-file exp_fun1.m, cmprt_ana_all.m and cmprt_plta.m
%             must be in the current directory or path.
%
%     07-Nov-2022 * Mack Gardner-Morse
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
resdir = fullfile('Results','All');    % Results directory
%
ifirst = true;          % First write to file
xlsnam = 'mri_fitra.xlsx';             % Results spreadsheet
xlsnam = fullfile(resdir,xlsnam);      % Include output directory
hdrs1 = {'Subject' 'Visit' 'Result' 'Leg' 'Load' 'Comprt' 'Bone' ...
         'Layer'};
hdrs2 = {'Pixels' 'T1R/T2S' 'RSS' 'ValidPix' 'Mean' 'Min' 'Max' ...
         'SD' 'COV'};
%
psnam = fullfile(resdir,'mri_fitra_'); % Start of PS file name
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
% Initialize Results Variables
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral, 2 = medial, and 3 = both
%   Index 6 - Bone - 1 = femur, 2 = tibia,  and 3 = both
%   Index 7 - Layer - 1 = deep, 2 = superficial, and both
%
t1r_res = zeros(nsubj,nvisit,2,2,3,3,3);
t1r_npx = zeros(nsubj,nvisit,2,2,3,3,3);
t1r_rss = zeros(nsubj,nvisit,2,2,3,3,3);
%
t1r_respx = cell(nsubj,nvisit,2,2,3,3,3);
t1r_rsspx = cell(nsubj,nvisit,2,2,3,3,3);
t1r_nps = cell(nsubj,nvisit,2,2,3,3,3);
%
t2s_res = zeros(nsubj,nvisit,2,2,3,3,3);
t2s_npx = zeros(nsubj,nvisit,2,2,3,3,3);
t2s_rss = zeros(nsubj,nvisit,2,2,3,3,3);
%
t2s_respx = cell(nsubj,nvisit,2,2,3,3,3);
t2s_rsspx = cell(nsubj,nvisit,2,2,3,3,3);
t2s_nps = cell(nsubj,nvisit,2,2,3,3,3);
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
      idr = contains(roinams,'2rois.mat','IgnoreCase',true);    % Small ROIs
      roinams = roinams(idr);          % ROI MAT files
      nroi = size(roinams,1);
%
      if nrho~=nroi
        error([' *** ERROR in mri_fitra:  Number of T1rho MAT', ...
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
         npix = prod(iszs);  % Number of pixels in an image
         fs = ['S' snt];     % Series number prefaced with a 'S'
%
         idm = contains(roinams,rhonam(1:end-4));     % Get matching file
         roinam = roinams{idm};
         load(fullfile(rdir,roinam),'icmprt','maskf','maskt','rsl');
%
% Parse Series Text for Leg and Load
%
         if strcmpi(st(1),'L')
           leg = 'L';
           legtxt = 'Left Leg ';
           ileg = 0;    % Coding for leg
         else
           leg = 'R';
           legtxt = 'Right Leg ';
           ileg = 1;
         end
%
         if contains(st,'Load','IgnoreCase',true)
           ld = 'LD';
           ldtxt = 'Loaded';
           ild = 1;     % Coding for load
         else
           ld = 'UL';
           ldtxt = 'Unloaded';
           ild = 0;
         end
%
         stitle = [legtxt ldtxt];      % Combine leg and load for plot titles
%
% Add Leg and Load to PS File Name
%
         psnamf = [psnamr leg '_' ld pstyp];     % Add leg and load to PS file name
%
% Combine Region Slices and Compartment Identifier
%
         idm = logical(icmprt-1);      % Index to medial compartment
%
         rsll = rsl(~idm);
         rslm = rsl(idm);
         rsls = {rsll; rslm; rsl};     % Lateral - row 1, medial - row 2, both - row 3
%          nrsls = cellfun('length',rsls);
%          nrsls = cellfun('size',rsls,1);
         nrsls = [size(rsll,1); size(rslm,1); size(rsl,1)];     % Faster than cellfun
%
% Combine Masks into a Cell Array
%
         maskf(:,2:3,:) = maskf;
         maskt(:,2:3,:) = maskt;
         maskf(:,1,:) = maskf(:,2,:)|maskf(:,3,:);    % Combine layers
         maskt(:,1,:) = maskt(:,2,:)|maskt(:,3,:);    % Combine layers
%
         falsel = false(npix,3,nrsls(1));   % Empty lateral mask
         falsem = false(npix,3,nrsls(2));   % Empty medial mask
%
         maskfl = maskf(:,:,~idm);     % Lateral
         maskfm = maskf(:,:,idm);      % Medial
%
         masktl = maskt(:,:,~idm);     % Lateral
         masktm = maskt(:,:,idm);      % Medial
%
         if icmprt(1)==1
           maskfl = cat(3,maskfl,falsem);
           masktl = cat(3,masktl,falsem);
%
           maskfm = cat(3,falsel,maskfm);
           masktm = cat(3,falsel,masktm);
         else
           maskfm = cat(3,maskfm,falsel);
           masktm = cat(3,masktm,falsel);
%
           maskfl = cat(3,falsem,maskfl);
           masktl = cat(3,falsem,masktl);
         end
%
         maskl = {maskfl; masktl; maskfl|masktl};     % Combine femur and tibia masks
         maskm = {maskfm; masktm; maskfm|masktm};     % Combine femur and tibia masks
         maskb = {maskf; maskt; maskf|maskt};         % Combine femur and tibia masks
%
         mask = {maskl; maskm; maskb};                % Combine compartment masks
%
% Do Compartmental Analysis
%
         [tc,~,rss,npx,id,tcp,ampp,rssp,nps] = cmprt_ana_all(v,mask, ...
                                 rsls,nrsls,splt,nslt,fun,init,tr0,opt);
         na = size(tc,1);              % Number of results
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral, 2 = medial, and 3 = both
%   Index 6 - Bone - 1 = femur, 2 = tibia, and 3 = both
%   Index 7 - Layer - 1 = deep, 2 = superficial, and 3 = both
%
% Note:  Layers for masks and compartment analysis variables are:
%        1 = both, 2 = superficial, and 3 = deep (opposite order).
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
         [~,idmx] = intersect(id,max(id),'row');
         sid = ['Subject ' subjnam];
         ttxt = {[sid ' on ' vnam]; stitle};
         cmprt_plta(v,maskb{3},rsls,nrsls,idt,tcp{idmx},nps{idmx}, ...
                    mxtr,cmap,ttxt,psnamf);
%
% Get Statistics on Pixel Results
%
         npxv = zeros(na,1);           % Number of valid results
         tcpm = zeros(na,1);           % Mean
         tcpmn = zeros(na,1);          % Minimum
         tcpmx = zeros(na,1);          % Maximum
         tcpsd = zeros(na,1);          % SD
%
         for ka = 1:na
            idv = tcp{ka}>=trmn&tcp{ka}<=trmx;
            npxv(ka) = sum(idv);       % Number of valid results
            tcpv = tcp{ka}(idv);       % Valid T1rho values
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
         ids = repmat(ids,na,1);
         ids = [ids id];               % All identifiers
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
      idr = contains(roinams,'2rois.mat','IgnoreCase',true);    % Small ROIs
      roinams = roinams(idr);          % ROI MAT files
      nroi = size(roinams,1);
%
      if nstar~=nroi
        error([' *** ERROR in mri_fitra:  Number of T2* MAT files', ...
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
         npix = prod(iszs);  % Number of pixels in an image
         fs = ['S' snt];     % Series number prefaced with a 'S'
%
         idm = contains(roinams,starnam(1:end-4));    % Get matching file
         roinam = roinams{idm};
         load(fullfile(rdir,roinam),'icmprt','maskf','maskt','rsl');
%
% Parse Series Text for Leg and Load
%
         if strcmpi(st(1),'L')
           leg = 'L';
           legtxt = 'Left Leg ';
           ileg = 0;    % Coding for leg
         else
           leg = 'R';
           legtxt = 'Right Leg ';
           ileg = 1;
         end
%
         if contains(st,'Load','IgnoreCase',true)
           ld = 'LD';
           ldtxt = 'Loaded';
           ild = 1;     % Coding for load
         else
           ld = 'UL';
           ldtxt = 'Unloaded';
           ild = 0;
         end
%
         stitle = [legtxt ldtxt];      % Combine leg and load for plot titles
%
% Add Leg and Load to PS File Name
%
         psnamf = [psnamr leg '_' ld pstyp];     % Add leg and load to PS file name
%
% Combine Region Slices and Compartment Identifier
%
         idm = logical(icmprt-1);      % Index to medial compartment
%
         rsll = rsl(~idm);
         rslm = rsl(idm);
         rsls = {rsll; rslm; rsl};     % Lateral - row 1, medial - row 2, both - row 3
%          nrsls = cellfun('length',rsls);
%          nrsls = cellfun('size',rsls,1);
         nrsls = [size(rsll,1); size(rslm,1); size(rsl,1)];     % Faster than cellfun
%
% Combine Masks into a Cell Array
%
         maskf(:,2:3,:) = maskf;
         maskt(:,2:3,:) = maskt;
         maskf(:,1,:) = maskf(:,2,:)|maskf(:,3,:);    % Combine layers
         maskt(:,1,:) = maskt(:,2,:)|maskt(:,3,:);    % Combine layers
%
         falsel = false(npix,3,nrsls(1));   % Empty lateral mask
         falsem = false(npix,3,nrsls(2));   % Empty medial mask
%
         maskfl = maskf(:,:,~idm);     % Lateral
         maskfm = maskf(:,:,idm);      % Medial
%
         masktl = maskt(:,:,~idm);     % Lateral
         masktm = maskt(:,:,idm);      % Medial
%
         if icmprt(1)==1
           maskfl = cat(3,maskfl,falsem);
           masktl = cat(3,masktl,falsem);
%
           maskfm = cat(3,falsel,maskfm);
           masktm = cat(3,falsel,masktm);
         else
           maskfm = cat(3,maskfm,falsel);
           masktm = cat(3,masktm,falsel);
%
           maskfl = cat(3,falsem,maskfl);
           masktl = cat(3,falsem,masktl);
         end
%
         maskl = {maskfl; masktl; maskfl|masktl};     % Combine femur and tibia masks
         maskm = {maskfm; masktm; maskfm|masktm};     % Combine femur and tibia masks
         maskb = {maskf; maskt; maskf|maskt};         % Combine femur and tibia masks
%
         mask = {maskl; maskm; maskb};                % Combine compartment masks
%
% Do Compartmental Analysis
%
         [tc,~,rss,npx,id,tcp,ampp,rssp,nps] = cmprt_ana_all(v,mask, ...
                                 rsls,nrsls,etns,netn,fun,init,ts0,opt);
         na = size(tc,1);              % Number of results
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral, 2 = medial, and 3 = both
%   Index 6 - Bone - 1 = femur, 2 = tibia, and 3 = both
%   Index 7 - Layer - 1 = deep, 2 = superficial, and 3 = both
%
% Note:  Layers for masks and compartment analysis variables are:
%        1 = both, 2 = superficial, and 3 = deep (opposite order).
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
         [~,idmx] = intersect(id,max(id),'row');
         sid = ['Subject ' subjnam];
         ttxt = {[sid ' on ' vnam]; stitle};
         cmprt_plta(v,maskb{3},rsls,nrsls,idt,tcp{idmx},nps{idmx}, ...
                    mxtr,cmap,ttxt,psnamf);
%
% Get Statistics on Pixel Results
%
         npxv = zeros(na,1);           % Number of valid results
         tcpm = zeros(na,1);           % Mean
         tcpmn = zeros(na,1);          % Minimum
         tcpmx = zeros(na,1);          % Maximum
         tcpsd = zeros(na,1);          % SD
%
         for ka = 1:na
            idv = tcp{ka}>=tsmn&tcp{ka}<=tsmx;
            npxv(ka) = sum(idv);       % Number of valid results
            tcpv = tcp{ka}(idv);       % Valid T2* values
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
         ids = repmat(ids,na,1);
         ids = [ids id];               % All identifiers
%
% Create and Write Table of Results
%
         t1 = array2table(ids,'VariableNames',hdrs1);
         t2 = table(npx,tc,rss,npxv,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
                    'VariableNames',hdrs2);
         t = [t1 t2];
%
         writetable(t,xlsnam,'WriteMode','append', ...
                    'WriteVariableNames',false);
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
save(fullfile(resdir,'mri_fitra.mat'),'t1r_res','t1r_npx','t1r_rss', ...
     't1r_respx','t1r_rsspx','t1r_nps','t2s_res','t2s_npx', ...
     't2s_rss','t2s_respx','t2s_rsspx','t2s_nps');
%
return