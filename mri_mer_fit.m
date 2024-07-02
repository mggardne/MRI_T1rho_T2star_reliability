%#######################################################################
%
%            * MRI Meniscus Erode Reliability FIT Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     MAT files and fits a monoexponential to the MRI data as a function
%     of spin lock or echo times.  T1rho or T2* are the time constants
%     of the fits.  The masks are eroded by one pixel on the boundary
%     using a cross, and square structural elements.  These masks are
%     used to get T1rho and T2* over the smaller areas.  Resulting
%     T1rho and T2* values and summary statistics are written to the
%     MS-Excel spreadsheets, mri_mer_fit_no.xlsx, mri_mer_fit_ec.xlsx,
%     and mri_mer_fit_es.xlsx, in the "Results" directory.
%
%     NOTES:  1.  Data MAT files must be in subdirectories "Visit1" or
%             "Visit2" in subject directories "MRIR *" where "*" is the
%             subject number.
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "_mrois".  See rd_dicom.m and
%             seg_mr_rois.m.
%
%             3.  M-file exp_fun1.m, cmprt_ana4mb.m, and cmprt_plt4mb.m
%             must be in the current directory or path.
%
%     26-Jun-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Setup Color Map
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
%
% Bounding Box for Plot Zoom Window
%
% boxr = [150 345 185 320];    % T1rho bounding box - not used - see bbox
% boxs = [90 250 150 250];     % T2* bounding box - not used - see bbox
%
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter', ...
               2e+3,'Algorithm','levenberg-marquardt','Jacobian', ...
               'on','UseParallel',true);
%
fun = @exp_fun1;        % Exponential function
%
% Get Different Structured Elements for Eroding the Masks by One Pixel
% (Increase 3 to 5 to erode by two (2) pixels.)
%
sec = strel(logical([0 1 0; 1 1 1; 0 1 0]));     % Crossed lines
ses = strel('square',3);                         % Square
%
ana_typ = ['no'; 'ec'; 'es'];          % Analysis identifiers
%
% Initialize Parameters
%
iskip = true;           % Skip T2* analysis
%
% init = -1;              % Use weighted least squares for starting parameters
% init = 0;               % Use linear least squares for starting parameters
init = 1;               % Use fixed starting parameters
tr0 = 20;               % 21.45 ms - mean w/ 50 ms threshold
% tr0 = 65;               % Initial T1rho estimate in ms
% tr0 = 80;               % Initial T1rho estimate in ms
trmx = 100;             % Maximum valid T1rho result
trmn = 0;               % Minimum valid T1rho result
ts0 = 12;               % 12.38 ms - mean w/ 50 ms threshold
% ts0 = 35;               % Initial T2* estimate in ms
tsmx = 100;             % Maximum valid T2* result
tsmn = 0;               % Minimum valid T2* result
%
mxtr = 60;              % Maximum scale on T1rho plots
mxts = 55;              % Maximum scale on T2* plots
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results','Meniscus');    % Results directory
xlsnam1 = 'mri_mer_fit_';              % Results spreadsheet
xlsnam1 = fullfile(resdir,xlsnam1);    % Include output directory
xlstyp = '.xlsx';       % Spreadsheet file extension
hdrs1 = {'Subject' 'Visit' 'Result' 'Leg' 'Load' 'Comprt' 'AP'};
hdrs2 = {'Pixels' 'T1R/T2S' 'RSS' 'ValidPix' 'Mean' 'Min' 'Max' ...
         'SD' 'COV'};
%
psnam = fullfile(resdir,'mri_mer_fit_');    % Start of PS file name
pstyp = '.ps';          % PS file type
%
% Get Subject Directories
%
spath = 'MRIR *';         % Path to series MAT files
sdirs = dir(spath);
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
nsubj = size(sdirs,1);
%
vdirs = {'Visit1'; 'Visit2'};
nvdir = 2;
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
%   Index 5 - Compartment - 1 = lateral, 2 = medial, and 3 = both
%   Index 6 - AP - 1 = anterior, 2 = posterior, and 3 = both
%
t1r_res = zeros(nsubj,2,2,2,3,3);
t1r_npx = zeros(nsubj,2,2,2,3,3);
t1r_rss = zeros(nsubj,2,2,2,3,3);
%
t1r_respx = cell(nsubj,2,2,2,3,3);
t1r_rsspx = cell(nsubj,2,2,2,3,3);
t1r_nps = cell(nsubj,2,2,2,3,3);
%
t2s_res = zeros(nsubj,2,2,2,3,3);
t2s_npx = zeros(nsubj,2,2,2,3,3);
t2s_rss = zeros(nsubj,2,2,2,3,3);
%
t2s_respx = cell(nsubj,2,2,2,3,3);
t2s_rsspx = cell(nsubj,2,2,2,3,3);
t2s_nps = cell(nsubj,2,2,2,3,3);
%
% Loop through Erosion Analysis with Different Masks
%
for ke = 1:3
%
% Add Erosion Analysis to File Names
%
   atyp = ana_typ(ke,:);
%
   psnama = [psnam atyp];              % Add analysis to PS file name
   xlsnam = [xlsnam1 atyp xlstyp];     % Add analysis to spreadsheet file name
%
   ifirst = true;       % First write to spreadsheet file
%
   atyptxt = [upper(atyp) ' Erosion']; % Text for plots
%
% Loop through Subjects
%
   for ks = 1:nsubj
%
% Get Subject Directory, Name and Number
%
      sdir = sdirs{ks};                % Current subject directory
      subj = eval(sdir(6:7));          % Subject number
      if subj<10
        snum = ['0' int2str(subj)];
      else
        snum = int2str(subj);
      end
%
      psnams = [psnama '_' snum];      % Add subject to PS file name
%
% Loop through Visits
%
      for kv = 1:nvdir
%
% Get Visit Subdirectory, Name and Number
%
         vdir = vdirs{kv};             % Current visit directory
         vstr = int2str(kv);           % Visit number as text
         vnam = ['Visit ' vstr];       % Visit name as text
         vid = kv-1;                   % Visit number
%
         psnamv = [psnams '_V' vstr];  % Add visit to PS file name
%
% Directory with Data Matlab MAT Files
%
         svdir = fullfile(sdir,vdir);
%
% Get T1rho MAT Files in Directory
%
         d = dir(fullfile(svdir,'T1rho_S*.mat'));
         rhonams = {d.name}';
         idr = contains(rhonams,'roi','IgnoreCase',true);  % ROI files
         rhonams = rhonams(~idr);
         idr = contains(rhonams,'chk','IgnoreCase',true);  % Check files
         rhonams = rhonams(~idr);
         id3 = contains(rhonams,'_3_');          % _3_?.mat files
         rhonams = rhonams(~id3);
         nrho = size(rhonams,1);
%
         d = dir(fullfile(svdir,'T1rho_S*_mrois.mat'));    % ROI files
         roinams = {d.name}';
         nroi = size(roinams,1);
%
         if nrho~=nroi
           error([' *** ERROR in mri_mer_fit:  Number of T1rho MAT', ...
                  ' files not equal to number of ROI MAT files!']);
         end
%
% T1rho Identifier
%
         ires = 0;      % ires = 0 - T1rho, ires = 1 - T2*
         idt = 1;       % Spin lock/echo time for plots - 1 = 0 ms spin lock time
%
         psnamr = [psnamv '_T1R_'];    % Add result type to PS file name
%
% Loop through T1rho MAT Files
%
         for km = 1:nroi
%
% Load Data
%
            roinam = roinams{km};
            load(fullfile(svdir,roinam),'bbox','cmprt','maskmal', ...
                 'maskmam','maskmpl','maskmpm','rsll','rslm');
%
            idm = contains(rhonams,roinam(1:end-10)); % Get matching file
            if ~any(idm)
              error([' *** ERROR in mri_mer_fit:  Matching T1rho', ...
                     ' MAT file not found for ROI MAT file:  ', ...
                     roinam '!']);
            end
            rhonam = rhonams{idm};
            load(fullfile(svdir,rhonam),'iszs','nslt','scmx','sns', ...
                 'snt','splt','st','v');
%
            npix = prod(iszs);    % Number of pixels in an image
            fs = ['S' snt];       % Series number prefaced with a 'S'
%
% Parse Series Text for Leg and Load
%
            if strcmpi(st(1),'L')
              leg = 'L';
              ileg = 0;    % Coding for left leg
              ltxt = 'Left';
            else
              leg = 'R';
              ileg = 1;    % Coding for right leg
              ltxt = 'Right';
            end
%
            if contains(st,'Load','IgnoreCase',true)
              ld = 'LD';               % Loaded
              iload = 1;               % Coding for loaded
              ldstr = 'Loaded';
            else
              ld = 'UL';               % Unloaded
              iload = 0;               % Coding for unloaded
              ldstr = 'Unloaded';
            end
%
% Add Leg and Load to PS File Name
%
            psnamf = [psnamr leg '_' ld pstyp];  % Add leg and load to PS file name
%
% Get Compartment Identifier for Slices
%
            rsllm = [rsll; rslm];      % Combined lateral and medial slices
            rsls = {rsll'; rslm'; rsllm'};  % Lateral - row 1, medial - row 2, both - row 3
            nrsls = [size(rsll,1); size(rslm,1); size(rsllm,1)];
%
% Get Masks
%
            emaskmal = maskmal;        % No erosion of the mask
            emaskmam = maskmam;        % No erosion of the mask
            emaskmpl = maskmpl;        % No erosion of the mask
            emaskmpm = maskmpm;        % No erosion of the mask
%
            if ke>1     % Erode masks
%
% Loop through the Lateral Slices
%
              for ksl = 1:nrsls(1)
%
% Get Lateral Slice Masks as Two-Dimensional Images
%
                 msk_tmpa = maskmal(:,ksl);
                 msk_tmpa = reshape(msk_tmpa,iszs);
                 msk_tmpp = maskmpl(:,ksl);
                 msk_tmpp = reshape(msk_tmpp,iszs);
%
% Erode Lateral Masks
%
                 if ke==2
                   msk_tmpa = imerode(msk_tmpa,sec);  % Crossed lines
                   msk_tmpp = imerode(msk_tmpp,sec);  % Crossed lines
                 else
                   msk_tmpa = imerode(msk_tmpa,ses);  % Square
                   msk_tmpp = imerode(msk_tmpp,ses);  % Square
                 end
                 emaskmal(:,ksl) =   msk_tmpa(:);
                 emaskmpl(:,ksl) =   msk_tmpp(:);
              end
%
% Loop through the Medial Slices
%
              for ksl = 1:nrsls(2)
%
% Get Medial Slice Masks as Two-Dimensional Images
%
                 msk_tmpa = maskmam(:,ksl);
                 msk_tmpa = reshape(msk_tmpa,iszs);
                 msk_tmpp = maskmpm(:,ksl);
                 msk_tmpp = reshape(msk_tmpp,iszs);
%
% Erode Medial Masks
%
                 if ke==2
                   msk_tmpa = imerode(msk_tmpa,sec);  % Crossed lines
                   msk_tmpp = imerode(msk_tmpp,sec);  % Crossed lines
                 else
                   msk_tmpa = imerode(msk_tmpa,ses);  % Square
                   msk_tmpp = imerode(msk_tmpp,ses);  % Square
                 end
                 emaskmam(:,ksl) =   msk_tmpa(:);
                 emaskmpm(:,ksl) =   msk_tmpp(:);
              end
%           
            end         % End of if - ke>1
%
% Combine Masks into a Cell Array
%
            maskmapl = {emaskmal; emaskmpl; emaskmal|emaskmpl}; % Combine lateral anterior and posterior masks
            maskmapm = {emaskmam; emaskmpm; emaskmam|emaskmpm}; % Combine medial anterior and posterior masks
            maskmalm = [emaskmal emaskmam]; % Combine lateral and medial anterior masks
            maskmplm = [emaskmpl emaskmpm]; % Combine lateral and medial posterior masks
            maskmaplm = [maskmapl{3} maskmapm{3}];    % Combine lateral and medial combined anterior and posterior masks
            maskmlm = {maskmalm; maskmplm; maskmaplm};     % Combine lateral and medial masks
            mask = {maskmapl; maskmapm; maskmlm};     % Combine compartment masks
%
% Do Compartmental Analysis
%
            [tc,~,rss,npx,id,tcp,ampp,rssp,nps] = cmprt_ana4mb(v, ...
                            mask,rsls,nrsls,splt,nslt,fun,init,tr0,opt);
            na = size(tc,1);           % Number of results
%
% Save Results
%
% Indices key:
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral, 2 = medial, and 3 = both
%   Index 6 - AP - 1 = anterior, 2 = posterior, and 3 = both
%
            for ka = 1:na
               t1r_res(ks,kv,ileg+1,iload+1,id(ka,1)+1,id(ka,2)+1) = ...
                                                                 tc(ka);
               t1r_npx(ks,kv,ileg+1,iload+1,id(ka,1)+1,id(ka,2)+1) ...
                                                              = npx(ka);
               t1r_rss(ks,kv,ileg+1,iload+1,id(ka,1)+1,id(ka,2)+1) = ...
                                                                rss(ka);
               t1r_respx{ks,kv,ileg+1,iload+1,id(ka,1)+1,id(ka,2)+1} ...
                                                              = tcp{ka};
               t1r_rsspx{ks,kv,ileg+1,iload+1,id(ka,1)+1,id(ka,2)+1} ...
                                                             = rssp{ka};
               t1r_nps{ks,kv,ileg+1,iload+1,id(ka,1)+1,id(ka,2)+1} = ...
                                                                nps{ka};
            end
%
% Plot Results
%
            sid = {['Subject ' snum ', ' vnam ', ' ltxt ' Leg']; ...
                   [ldstr ', T1\rho, ' atyptxt]};
            cmprt_plt4mb(v,mask,rsls,nrsls,idt,tcp,nps,mxtr,cmap, ...
                         bbox,sid,psnamf);
%
% Get Statistics on Pixel Results
%
            npxv = zeros(na,1);        % Number of valid results
            tcpm = zeros(na,1);        % Mean
            tcpmn = zeros(na,1);       % Minimum
            tcpmx = zeros(na,1);       % Maximum
            tcpsd = zeros(na,1);       % SD
%
            for ka = 1:na
               idv = tcp{ka}>=trmn&tcp{ka}<=trmx;
               npxv(ka) = sum(idv);         % Number of valid results
               if npxv(ka)>0                
                 tcpv = tcp{ka}(idv);       % Valid T1rho values
                 tcpm(ka) = mean(tcpv);     % Mean
                 tcpmn(ka) = min(tcpv);     % Minimum
                 tcpmx(ka) = max(tcpv);     % Maximum
                 tcpsd(ka) = std(tcpv);     % SD
               end
            end
%
            tcpcov = 100*tcpsd./tcpm;  % Coefficient of variation
            tcpcov(isnan(tcpcov)) = 0; % Catch any NaNs
%
% Combine Identifiers
%
            ids = [subj vid ires ileg iload];    % MAT file identifiers
            ids = repmat(ids,na,1);
            ida = [ids id];            % All identifiers
%
% Create and Write Table of Results
%
            t1 = array2table(ida,'VariableNames',hdrs1);
            t2 = table(npx,tc,rss,npxv,tcpm,tcpmn,tcpmx,tcpsd, ...
                       tcpcov,'VariableNames',hdrs2);
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
         end            % End of km loop - T1rho MAT file loop
%
         close all;     % Close all plot windows
%
% Get T2* MAT Files in Directory
%
         if ~iskip
      d = dir(fullfile(spdir,'T2star_S*.mat'));
      starnams = {d.name}';
      idr = contains(starnams,'roi','IgnoreCase',true);    % ROI files
      starnams = starnams(~idr);
      idr = contains(starnams,'chk','IgnoreCase',true);    % Check files
      starnams = starnams(~idr);
      nstar = size(starnams,1);
%
      d = dir(fullfile(sdir,'T2star_S*mrois.mat'));
      roinams = {d.name}';
      nroi = size(roinams,1);
%
      if nstar~=nroi
        error([' *** ERROR in mri_mer_fit:  Number of T2* MAT', ...
               ' files not equal to number of ROI MAT files!']);
      end
%
% T2* Identifier
%
      ires = 1;         % ires = 0 - T1rho, ires = 1 - T2*
      idt = 3;          % Spin lock/echo time for plots - 3 = 5 ms echo time
%
      psnamr = [psnams '_T2S_'];       % Add result type to PS file name
%
% Loop through T2* MAT Files
%
      for km = 1:nroi
%
% Load Data
%
         roinam = roinams{km};
         load(fullfile(sdir,roinam),'bbox','cmprt','maskmal', ...
              'maskmam','maskmpl','maskmpm','rsll','rslm');
%
         idm = contains(starnams,roinam(1:end-10));   % Get matching file
         if ~any(idm)
           error([' *** ERROR in mri_mer_fit:  Matching T2* MAT', ...
                  ' file not found for ROI MAT file:  ' roinam '!']);
         end
         starnam = starnams{idm};
         load(fullfile(spdir,starnam),'etns','iszs','netn','scmx', ...
              'sns','snt','st','v');
%
         npix = prod(iszs);  % Number of pixels in an image
         fs = ['S' snt];     % Series number prefaced with a 'S'
%
% Parse Series Text for Leg
%
         if strcmpi(st(1),'L')
           leg = 'L';
           ileg = 0;       % Coding for left leg
           ltxt = 'Left';
         else
           leg = 'R';
           ileg = 1;       % Coding for right leg
           ltxt = 'Right';
         end
%
% Add Leg to PS File Name
%
         psnamf = [psnamr leg pstyp];  % Add leg to PS file name
%
% Get Compartment Identifier for Slices
%
         rsls = {rsll'; rslm'};        % Lateral - row 1, medial - row 2
         nrsls = [size(rsll,1); size(rslm,1)];
%
% Get Masks
%
         emaskmal = maskmal;           % No erosion of the mask
         emaskmam = maskmam;           % No erosion of the mask
         emaskmpl = maskmpl;           % No erosion of the mask
         emaskmpm = maskmpm;           % No erosion of the mask
%
         if ke>1       % Erode masks
%
% Loop through the Lateral Slices
%
           for ksl = 1:nrsls(1)
%
% Get Lateral Slice Masks as Two-Dimensional Images
%
              msk_tmpa = maskmal(:,ksl);
              msk_tmpa = reshape(msk_tmpa,iszs);
              msk_tmpp = maskmpl(:,ksl);
              msk_tmpp = reshape(msk_tmpp,iszs);
%
% Erode Lateral Masks
%
              if ke==2
                msk_tmpa = imerode(msk_tmpa,sec);  % Crossed lines
                msk_tmpp = imerode(msk_tmpp,sec);  % Crossed lines
              else
                msk_tmpa = imerode(msk_tmpa,ses);  % Square
                msk_tmpp = imerode(msk_tmpp,ses);  % Square
              end
              emaskmal(:,ksl) =   msk_tmpa(:);
              emaskmpl(:,ksl) =   msk_tmpp(:);
           end
%
% Loop through the Medial Slices
%
           for ksl = 1:nrsls(2)
%
% Get Medial Slice Masks as Two-Dimensional Images
%
              msk_tmpa = maskmam(:,ksl);
              msk_tmpa = reshape(msk_tmpa,iszs);
              msk_tmpp = maskmpm(:,ksl);
              msk_tmpp = reshape(msk_tmpp,iszs);
%
% Erode Medial Masks
%
              if ke==2
                msk_tmpa = imerode(msk_tmpa,sec);  % Crossed lines
                msk_tmpp = imerode(msk_tmpp,sec);  % Crossed lines
              else
                msk_tmpa = imerode(msk_tmpa,ses);  % Square
                msk_tmpp = imerode(msk_tmpp,ses);  % Square
              end
              emaskmam(:,ksl) =   msk_tmpa(:);
              emaskmpm(:,ksl) =   msk_tmpp(:);
           end
%
         end            % End of if - ke>1
%
% Combine Masks into a Cell Array
%
         maskmapl = {emaskmal; emaskmpl};   % Combine lateral anterior and posterior masks
         maskmapm = {emaskmam; emaskmpm};   % Combine medial anterior and posterior masks
         mask = {maskmapl; maskmapm};       % Combine compartment masks
%
% Do Compartmental Analysis
%
         [tc,~,rss,npx,id,tcp,ampp,rssp,nps] = cmprt_ana4m(v,mask, ...
                                    rsls,nrsls,etns,netn,fun,init,ts0,opt);
         na = size(tc,1);              % Number of results
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Compartment - 1 = lateral and 2 = medial
%   Index 4 - AP - 1 = anterior and 2 = posterior
%
         for ka = 1:na
            t2s_res(ks,ileg+1,id(ka,1)+1,id(ka,2)+1) = tc(ka);
            t2s_npx(ks,ileg+1,id(ka,1)+1,id(ka,2)+1) = npx(ka);
            t2s_rss(ks,ileg+1,id(ka,1)+1,id(ka,2)+1) = rss(ka);
            t2s_respx{ks,ileg+1,id(ka,1)+1,id(ka,2)+1} = tcp{ka};
            t2s_rsspx{ks,ileg+1,id(ka,1)+1,id(ka,2)+1} = rssp{ka};
            t2s_nps{ks,ileg+1,id(ka,1)+1,id(ka,2)+1} = nps{ka};
         end
%
% Plot Results
%
         sid = ['Subject ' sdir ', ' ltxt ' Leg, T2*, ' atyptxt];
         cmprt_plt4m(v,mask,rsls,nrsls,idt,tcp,nps,mxts,cmap,bbox, ...
                     sid,psnamf);
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
            if npxv(ka)>0
              tcpv = tcp{ka}(idv);     % Valid T2* values
              tcpm(ka) = mean(tcpv);   % Mean
              tcpmn(ka) = min(tcpv);   % Minimum
              tcpmx(ka) = max(tcpv);   % Maximum
              tcpsd(ka) = std(tcpv);   % SD
            end
         end
%
         tcpcov = 100*tcpsd./tcpm;     % Coefficient of variation
         tcpcov(isnan(tcpcov)) = 0;    % Catch any NaNs
%
% Combine Identifiers
%
         ids = [subj ires ileg];       % MAT file identifiers
         ids = repmat(ids,na,1);
         ida = [ids id];               % All identifiers
%
% Create and Write Table of Results
%
         t1 = array2table(ida,'VariableNames',hdrs1);
         t2 = table(npx,tc,rss,npxv,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
                    'VariableNames',hdrs2);
         t = [t1 t2];
%
         writetable(t,xlsnam,'WriteMode','append', ...
                    'WriteVariableNames',false);
%
      end               % End of km loop - T2* MAT file loop
%
           close all;   % Close all plot windows
%
         end            % End of if - iskip
%
      end               % End of kv loop - visits loop
%
   end                  % End of ks loop - subjects loop
%
% Save to MAT File
%
   matnam = fullfile(resdir,['mri_mer_fit_' atyp '.mat']);
   if ~iskip
     save(matnam,'t1r_res','t1r_npx','t1r_rss','t1r_respx', ...
          't1r_rsspx','t1r_nps','t2s_res','t2s_npx','t2s_rss', ...
          't2s_respx','t2s_rsspx','t2s_nps');
   else
     save(matnam,'t1r_res','t1r_npx','t1r_rss','t1r_respx', ...
          't1r_rsspx','t1r_nps');
   end
%
end                     % End of ke loop - erosion analysis loop
%
return