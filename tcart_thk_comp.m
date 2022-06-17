%##########################################################################
%
%       *Tibial Cartilage Thickness Comparison Program*
%
%           M-file used to compare two tibial cartilage thickness maps
%           generated and stored in the MRIR database
%      
%
%   Notes:
%           1. Will take you direclty to MRIR directory to find files
%
%           2. Order to selecting files in selection Window: 
%              Subject Number> Visit 2 > Scan (FFE vs RHO)  
%
%           3.
%
%
%
%
%
%
%           16-March-2022  *  Andy "Big Brain" Borah
%
%
%##########################################################################

%% Plotting on 
iplt = true;  % on
%iplt = false;  % off

% plot printing
iprt = true;  % on
%iprt= false:  % off

%% Directory Selection

ddir_ffe = uigetdir('L:\Labs\Beynnon\General\Regulatory\Protocols\STUDY00001137_BeynnonFiorentinoTourville_MRIReliability\Working Binder\SUBJECT DATA','Select FFE Scan To Analyze');
fold_ffe = dir(ddir_ffe);
mat_ffe = fullfile(ddir_ffe,'tcart08_1.mat');

ddir_rho = uigetdir('L:\Labs\Beynnon\General\Regulatory\Protocols\STUDY00001137_BeynnonFiorentinoTourville_MRIReliability\Working Binder\SUBJECT DATA','Select RHO Scan To Analyze');
fold_rho = dir(ddir_rho);
mat_rho = fullfile(ddir_rho,'tcart08_1.mat');
%% Load in Variables and Assign Unique Names

load(mat_ffe)
cthkls_ffe = cthkls;
cthkms_ffe = cthkms;
ilegs_ffe = ilegs;
kids_ffe = kids;
snams_ffe = snams;
tribl_ffe = tribl;          % ffe Variables
tribm_ffe = tribm;
trilss_ffe = trilss;
trimss_ffe = trimss;
xyzbl_ffe = xyzbl;
xyzbm_ffe = xyzbm;
xyzils_ffe = xyzils;
xyzims_ffe = xyzims;
xyzlss_ffe = xyzlss;
xyzmss_ffe = xyzmss;
zgl_ffe = zgl;
zgm_ffe = zgm;

load(mat_rho)
cthkls_rho = cthkls;
cthkms_rho = cthkms;
ilegs_rho = ilegs;
kids_rho = kids;
snams_rho = snams;
tribl_rho = tribl;          % Rho Variables
tribm_rho = tribm;
trilss_rho = trilss;
trimss_rho = trimss;
xyzbl_rho = xyzbl;
xyzbm_rho = xyzbm;
xyzils_rho = xyzils;
xyzims_rho = xyzims;
xyzlss_rho = xyzlss;
xyzmss_rho = xyzmss;
zgl_rho = zgl;
zgm_rho = zgm;

clear cthkls;
clear cthkms;
clear ilegs;
clear kids;
clear snams;
clear tribl;            % Clearing Workspace of Unused Variables
clear tribm;
clear trilss;
clear trimss;
clear xyzbl;
clear xyzbm;
clear xyzils;
clear xyzims;
clear xyzlss;
clear xyzmss;
clear zgl;
clear zgm;
%% Plot Raw Data

%% Calculate Differences

% Lateral Compartment
% Add Zeroes to Create Same Size Arrays
sffe = size(cthkls_ffe);
srho = size(cthkls_rho);
sdiff = sffe(1)-srho(1);
if sdiff > 0
    cthkls_rho(srho(1)+sdiff,1) = 0;  % Adding Zeroes to the Rho Array
end
if sdiff < 0
    cthkls_ffe(sffe(1)+abs(sdiff)) = 0;  % Adding Zeroes to FFE Array
end

sffe = size(cthkls_ffe);

numel = sffe(1);
lat_diff_calc = zeros(numel,1);
for k = 1:numel
    if isnan(cthkls_ffe(k,1))
        cthkls_ffe(k,1) = 0;
    end
    if isnan(cthkls_rho(k,1))
        cthkls_rho(k,1) = 0;
    end
    lat_diff_calc(k,1) = (cthkls_ffe(k,1)-cthkls_rho(k,1));
end

% Medial Compartment
% Add Zeroes to Create Same Size Arrays
sffe = size(cthkms_ffe);
srho = size(cthkms_rho);
sdiff = sffe(1)-srho(1);
if sdiff > 0
    cthkms_rho(srho(1)+sdiff,1) = 0;  % Adding Zeroes to the Rho Array
end
if sdiff < 0
    cthkms_ffe(sffe(1)+abs(sdiff)) = 0;  % Adding Zeroes to FFE Array
end

sffe = size(cthkms_ffe);

numel = sffe(1);
med_diff_calc = zeros(numel,1);

for k = 1:numel
    if isnan(cthkms_ffe(k,1))
        cthkms_ffe(k,1) = 0;
    end
    if isnan(cthkms_rho(k,1))
        cthkms_rho(k,1) = 0;
    end
    med_diff_calc(k,1) = (cthkms_ffe(k,1)-cthkms_rho(k,1));
end
%% Plot

%% Save

rnam = fullfile(ddir_ffe,['tcart_thk_comp' '.mat']);    % Results MAT file
save(rnam,'med_diff_calc','lat_diff_calc');