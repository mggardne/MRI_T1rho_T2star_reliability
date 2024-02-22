%#######################################################################
%
%            * MRI FIT Reliability 4 Unloaded PLoT Program *
%
%          M-File which reads the registered MRI data, segmentation,
%     and analysis results MAT files and plots the T1rho and T2*
%     results for particular slices for subject 01-BC on visit 2.  The
%     resulting T1rho and T2* plots are written to JPEG files:
%     mri_fitr4u_t1r.jpg, and mri_fitr4u_t2s.jpg.
%
%     NOTES:  1.  Data MAT files must be in subject directory
%             "MRIR 01-BC" and visit subdirectory "Visit2".
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "2rois".  See rd_dicom.m and
%             seg_2rois.m.
%
%             3.  M-file cmprt_plt4us.m must be in the current 
%             directory or path.
%
%             4.  The analysis results MAT file, mri_fitr4u.mat, must
%             be in the subdirectory, "Results\Unloaded_Femur4" under
%             the main MRI_Reliability_Study\ folder.
%
%             5.  This program must start from the directory
%             "MRIR 01-BC\Visit2".
%
%             6.  Slices are identified by variable, "slk".
%
%             7.  The following letter(s) in variable names have the
%             following meanings:
%
%             [1,2] - 1 is region 1 - loaded contact region and 2 is
%             region 2 - unloaded posterior femur region.
%
%             [l,m] - l is the lateral compartment and m is the medial
%             compartment.
%
%             [f,t] - f is the femur and t is the tibia.
%
%             [r,s] - r is T1rho and s is T2*.
%             [t1r,t2s] - t1r is T1rho and t2s is T2*.
%
%             [t1r,t1ru] - no identifier is region 1 and u is region 2.
%             [t2s,t2su] - no identifier is region 1 and u is region 2.
%
%             8.  MRI series for 01-BC Left Leg:
%             T1rho Unloaded = Series 0301
%             T1rho Loaded = Series 1001
%             T2* Unloaded = Series 0401
%             T2* Loaded = Series 1201
%
%             9.  T1rho slice 20 is in the lateral compartment.  T2*
%             slices 15 and 16 are in the lateral compartment.
%
%     07-Dec-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Setup Color Map
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
%
% Initialize Parameters
%
trmx = 100;             % Maximum valid T1rho result
trmn = 0;               % Minimum valid T1rho result
tsmx = 100;             % Maximum valid T2* result
tsmn = 0;               % Minimum valid T2* result
%
mxtr = 80;              % Maximum scale on T1rho plots
mxts = 75;              % Maximum scale on T2* plots
%
% Output Directory and JPEG Output Files
%
resdir = fullfile('..\..\','Results','Unloaded_Femur4');   % Results directory
%
psnam = fullfile(resdir,'mri_fitr4u_');     % Start of JPEG file names
jext = '.jpg';          % JPEG file type
%
% Load Results
%
%       t1r_nps: {6-D cell}
%       t1r_npx: [6-D double]
%       t1r_res: [6-D double]
%     t1r_respx: {6-D cell}
%       t1r_rss: [6-D double]
%     t1r_rsspx: {6-D cell}
%      t1ru_nps: {5-D cell}
%      t1ru_npx: [5-D double]
%      t1ru_res: [5-D double]
%    t1ru_respx: {5-D cell}
%      t1ru_rss: [5-D double]
%    t1ru_rsspx: {5-D cell}
%       t2s_nps: {6-D cell}
%       t2s_npx: [6-D double]
%       t2s_res: [6-D double]
%     t2s_respx: {6-D cell}
%       t2s_rss: [6-D double]
%     t2s_rsspx: {6-D cell}
%      t2su_nps: {5-D cell}
%      t2su_npx: [5-D double]
%      t2su_res: [5-D double]
%    t2su_respx: {5-D cell}
%      t2su_rss: [5-D double]
%    t2su_rsspx: {5-D cell}
%
load(fullfile(resdir,'mri_fitr4u.mat'),'t1r_nps','t1r_respx', ...
     't1ru_nps','t1ru_respx','t2s_nps','t2s_respx', ...
     't2su_nps','t2su_respx');
%
% Unloaded T1rho Data
%
idt = 1;                % Spin lock/echo time for plots - 1 = 0 ms spin lock time
%
% Load Image Data
%
%     ddir: 's0301'
%    fnams: {256×1 cell}
%     iszs: [512 512]
%    nfile: 256
%     nsls: 64
%     nslt: 4
%    pspcs: [0.4883 0.4883]
%     scmx: 2540
%      sns: 301
%      snt: '301'
%     splt: [4×1 double]
%       st: 'L T1rhoL 0 10 40 80ms'
%        v: [512×512×64×4 double]
%
load('T1rho_S301.mat','v');
%
% Load Masks
%
%       brois: [1×2 struct]
%           f: {2×34 cell}
%          f3: {2×34 cell}
%       ibone: [34×2 logical]
%      icmprt: [34×1 double]
%       maskf: [262144×2×34 logical]
%    maskfrl1: [262144×2×8 logical]
%    maskfrl2: [262144×2×8 logical]
%    maskfrm1: [262144×2×8 logical]
%    maskfrm2: [262144×2×8 logical]
%      maskr1: [262144×2×34 logical]
%      maskr2: [262144×2×34 logical]
%       maskt: [262144×2×34 logical]
%    masktrl1: [262144×2×8 logical]
%    masktrm1: [262144×2×8 logical]
%         rsl: [34×1 double]
%       rsll1: [8×1 double]
%       rsll2: [8×1 double]
%       rslm1: [8×1 double]
%       rslm2: [8×1 double]
%           t: {2×34 cell}
%          t3: {2×34 cell}
%
load('T1rho_S301_2rois.mat','maskfrl1','maskfrl2','maskfrm1', ...
     'maskfrm2','masktrl1','masktrm1','rsll1','rsll2','rslm1','rslm2');
%
% Combine Masks into a Cell Array
%
maskl1 = {maskfrl1; masktrl1};         % Combine femur and tibia masks
maskm1 = {maskfrm1; masktrm1};         % Combine femur and tibia masks
mask1 = {maskl1; maskm1};              % Combine compartment masks
%
mask2 = {maskfrl2; maskfrm2};          % Combine compartment masks
%
% Combine Compartment Identifiers into Region Slices
%
rsls1 = {rsll1; rslm1};                % Lateral - row 1, medial - row 2
rsls2 = {rsll2; rslm2};                % Lateral - row 1, medial - row 2
%
% Get Unloaded T1rho Results
%
% Region 1 (contact femur and tibia cartilage)
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - Bone - 1 = femur and 2 = tibia
%
nps1 = squeeze(t1r_nps(1,2,1,1,1,:));
tcp1 = squeeze(t1r_respx(1,2,1,1,1,:));
%
% Region 2 (non-contact femur cartilage)
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%
nps2 = squeeze(t1ru_nps(1,2,1,1,:));
tcp2 = squeeze(t1ru_respx(1,2,1,1,:));
%
% Unloaded T1rho Image (Image a)
%
slk = 20;
rimgra = cmprt_plt4us(v,mask1,rsls1,mask2,rsls2,idt,slk,tcp1,nps1, ...
                      tcp2,nps2,mxtr);
%
% Loaded T1rho Data
%
% Load Image Data
%
%     ddir: 's1001'
%    fnams: {256×1 cell}
%     iszs: [512 512]
%    nfile: 256
%     nsls: 64
%     nslt: 4
%    pspcs: [0.4883 0.4883]
%     scmx: 2540
%      sns: 1001
%      snt: '1001'
%     splt: [4×1 double]
%       st: 'L T1rho Load_TSL 0 10 40 80ms'
%        v: [512×512×64×4 double]
%
load('T1rho_S1001.mat','v');
%
% Load Masks
%
%       brois: [1×2 struct]
%           f: {2×31 cell}
%          f3: {2×31 cell}
%       ibone: [31×2 logical]
%      icmprt: [31×1 double]
%       maskf: [262144×2×31 logical]
%    maskfrl1: [262144×2×8 logical]
%    maskfrl2: [262144×2×9 logical]
%    maskfrm1: [262144×2×8 logical]
%    maskfrm2: [262144×2×9 logical]
%      maskr1: [262144×2×31 logical]
%      maskr2: [262144×2×31 logical]
%       maskt: [262144×2×31 logical]
%    masktrl1: [262144×2×8 logical]
%    masktrm1: [262144×2×8 logical]
%         rsl: [31×1 double]
%       rsll1: [8×1 double]
%       rsll2: [9×1 double]
%       rslm1: [8×1 double]
%       rslm2: [9×1 double]
%           t: {2×31 cell}
%          t3: {2×31 cell}
%
load('T1rho_S1001_2rois.mat','maskfrl1','maskfrl2','maskfrm1', ...
     'maskfrm2','masktrl1','masktrm1','rsll1','rsll2','rslm1','rslm2');
%
% Combine Masks into a Cell Array
%
maskl1 = {maskfrl1; masktrl1};         % Combine femur and tibia masks
maskm1 = {maskfrm1; masktrm1};         % Combine femur and tibia masks
mask1 = {maskl1; maskm1};              % Combine compartment masks
%
mask2 = {maskfrl2; maskfrm2};          % Combine compartment masks
%
% Combine Compartment Identifiers into Region Slices
%
rsls1 = {rsll1; rslm1};                % Lateral - row 1, medial - row 2
rsls2 = {rsll2; rslm2};                % Lateral - row 1, medial - row 2
%
% Get Loaded T1rho Results
%
% Region 1 (contact femur and tibia cartilage)
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - Bone - 1 = femur and 2 = tibia
%
nps1 = squeeze(t1r_nps(1,2,1,2,1,:));
tcp1 = squeeze(t1r_respx(1,2,1,2,1,:));
%
% Region 2 (non-contact femur cartilage)
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%
nps2 = squeeze(t1ru_nps(1,2,1,2,:));
tcp2 = squeeze(t1ru_respx(1,2,1,2,:));
%
% Loaded T1rho Image (Image b)
%
rimgrb = cmprt_plt4us(v,mask1,rsls1,mask2,rsls2,idt,slk,tcp1,nps1, ...
                      tcp2,nps2,mxtr);
%
% Crop and Plot T1rho Unloaded and Loaded
%
xcropa = 120:310;       % Image height
ycropa = 135:350;       % Image width
%
xcropb = 126:316;       % Image height
ycropb = 138:353;       % Image width
%
rimgr = [rimgra(xcropa,ycropa) rimgrb(xcropb,ycropb)];
%
figure;
orient landscape;
%
imshow(rimgr,[-mxtr mxtr]);
colormap(cmap);
hc = colorbar;
hc.Limits = [0 mxtr];
%
print('-djpeg95','-r600',[psnam 't1r' jext]);
%
% Unloaded T2* Data
%
idt = 3;                % Spin lock/echo time for plots - 3 = 5 ms echo time
%
% Load Image Data
%
%     ddir: {5×1 cell} {'s0401'} {'s0501'} {'s0601'} {'s0701'} {'s0801'}
%     etns: [5×1 double] [0.4250 1.0000 5.0000 15.0010 30.0000]
%    fnams: {48×1 cell}
%      id5: [4×1 double] [3 3 3 3]
%     iszs: [320 320]
%     netn: 5
%    nfile: 240
%     nsls: 48
%    pspcs: [0.4375 0.4375]
%     scmx: 38580
%      sns: [5×1 double] [401 501 601 701 801]
%      snt: '401'
%       st: 'L TE0.42_T2* MAP'
%        v: [320×320×48×5 double]
%
load('T2star_S401.mat','v');
%
% Load Masks
%
%       brois: [1×2 struct]
%           f: {2×21 cell}
%          f3: {2×21 cell}
%       ibone: [21×2 logical]
%      icmprt: [21×1 double]
%       maskf: [102400×2×21 logical]
%    maskfrl1: [102400×2×6 logical]
%    maskfrl2: [102400×2×6 logical]
%    maskfrm1: [102400×2×6 logical]
%    maskfrm2: [102400×2×6 logical]
%      maskr1: [102400×2×21 logical]
%      maskr2: [102400×2×21 logical]
%       maskt: [102400×2×21 logical]
%    masktrl1: [102400×2×6 logical]
%    masktrm1: [102400×2×6 logical]
%         rsl: [21×1 double]
%       rsll1: [6×1 double]
%       rsll2: [6×1 double]
%       rslm1: [6×1 double]
%       rslm2: [6×1 double]
%           t: {2×21 cell}
%          t3: {2×21 cell}
%
load('T2star_S401_2rois.mat','maskfrl1','maskfrl2','maskfrm1', ...
     'maskfrm2','masktrl1','masktrm1','rsll1','rsll2','rslm1','rslm2');
%
% Combine Masks into a Cell Array
%
maskl1 = {maskfrl1; masktrl1};         % Combine femur and tibia masks
maskm1 = {maskfrm1; masktrm1};         % Combine femur and tibia masks
mask1 = {maskl1; maskm1};              % Combine compartment masks
%
mask2 = {maskfrl2; maskfrm2};          % Combine compartment masks
%
% Combine Compartment Identifiers into Region Slices
%
rsls1 = {rsll1; rslm1};                % Lateral - row 1, medial - row 2
rsls2 = {rsll2; rslm2};                % Lateral - row 1, medial - row 2
%
% Get Unloaded T2* Results
%
% Region 1 (contact femur and tibia cartilage)
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - Bone - 1 = femur and 2 = tibia
%
nps1 = squeeze(t2s_nps(1,2,1,1,1,:));
tcp1 = squeeze(t2s_respx(1,2,1,1,1,:));
%
% Region 2 (non-contact femur cartilage)
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%
nps2 = squeeze(t2su_nps(1,2,1,1,:));
tcp2 = squeeze(t2su_respx(1,2,1,1,:));
%
% Unloaded T2* Image (Image c)
%
slk = 15;
rimgsc = cmprt_plt4us(v,mask1,rsls1,mask2,rsls2,idt,slk,tcp1,nps1, ...
                      tcp2,nps2,mxts);
%
% Loaded T2* Data
%
% Load Image Data
%
%     ddir: {5×1 cell} {'s1201'} {'s1301'} {'s1401'} {'s1501'} {'s1601'}
%     etns: [5×1 double] [0.4250 1.0000 5.0000 15.0010 30.0000]
%    fnams: {48×1 cell}
%      id5: [4×1 double] [3 3 3 3]
%     iszs: [320 320]
%     netn: 5
%    nfile: 240
%     nsls: 48
%    pspcs: [0.4375 0.4375]
%     scmx: 40540
%      sns: [5×1 double] [1201 1301 1401 1501 1601]
%      snt: '1201'
%       st: 'L Load TE0.42_T2*map'
%        v: [320×320×48×5 double]
%
load('T2star_S1201.mat','v');
%
% Load Masks
%
%       brois: [1×2 struct]
%           f: {2×20 cell}
%          f3: {2×20 cell}
%       ibone: [20×2 logical]
%      icmprt: [20×1 double]
%       maskf: [102400×2×20 logical]
%    maskfrl1: [102400×2×6 logical]
%    maskfrl2: [102400×2×6 logical]
%    maskfrm1: [102400×2×6 logical]
%    maskfrm2: [102400×2×6 logical]
%      maskr1: [102400×2×20 logical]
%      maskr2: [102400×2×20 logical]
%       maskt: [102400×2×20 logical]
%    masktrl1: [102400×2×6 logical]
%    masktrm1: [102400×2×6 logical]
%         rsl: [20×1 double]
%       rsll1: [6×1 double]
%       rsll2: [6×1 double]
%       rslm1: [6×1 double]
%       rslm2: [6×1 double]
%           t: {2×20 cell}
%          t3: {2×20 cell}
%
load('T2star_S1201_2rois.mat','maskfrl1','maskfrl2','maskfrm1', ...
     'maskfrm2','masktrl1','masktrm1','rsll1','rsll2','rslm1','rslm2');
%
% Combine Masks into a Cell Array
%
maskl1 = {maskfrl1; masktrl1};         % Combine femur and tibia masks
maskm1 = {maskfrm1; masktrm1};         % Combine femur and tibia masks
mask1 = {maskl1; maskm1};              % Combine compartment masks
%
mask2 = {maskfrl2; maskfrm2};          % Combine compartment masks
%
% Combine Compartment Identifiers into Region Slices
%
rsls1 = {rsll1; rslm1};                % Lateral - row 1, medial - row 2
rsls2 = {rsll2; rslm2};                % Lateral - row 1, medial - row 2
%
% Get Loaded T2* Results
%
% Region 1 (contact femur and tibia cartilage)
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - Bone - 1 = femur and 2 = tibia
%
nps1 = squeeze(t2s_nps(1,2,1,2,1,:));
tcp1 = squeeze(t2s_respx(1,2,1,2,1,:));
%
% Region 2 (non-contact femur cartilage)
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Load - 1 = unloaded and 2 = loaded
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%
nps2 = squeeze(t2su_nps(1,2,1,2,:));
tcp2 = squeeze(t2su_respx(1,2,1,2,:));
%
% Loaded T2* Image (Image d)
%
slk = 16;
rimgsd = cmprt_plt4us(v,mask1,rsls1,mask2,rsls2,idt,slk,tcp1,nps1, ...
                      tcp2,nps2,mxts);
%
% Crop and Plot T1rho Unloaded and Loaded
%
xcropc = 50:255;        % Image height
ycropc = 40:275;        % Image width
%
xcropd = 57:262;        % Image height
ycropd = 39:274;        % Image width
%
rimgs = [rimgsc(xcropc,ycropc) rimgsd(xcropd,ycropd)];
%
figure;
orient landscape;
%
imshow(rimgs,[-mxts mxts]);
colormap(cmap);
hc = colorbar;
hc.Limits = [0 mxts];
%
print('-djpeg95','-r600',[psnam 't2s' jext]);
%
return