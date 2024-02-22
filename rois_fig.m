%
% Radius of Analysis Region and Segmentation Parameter
%
rrad = 6.0;             % Radius of analysis region (width of 3 or 4 slices)
rrad2 = rrad*rrad;      % Radius squared
itroch = true;          % Read trochlea?
% itroch = false;         % Read trochlea?
thresh = 1.953125;      % Four (4) pixels = 4*0.48828125 = 1.953125
% iplt = true;            % Plot slices
iplt = false;           % No slice plots
% vw = [127.5 30];        % Original view
vw = [120 37.5];        % Better view of contact ROIs?
% vw = [140 40];        % Better view of contact ROIs?
%
% Axis Vector and Contact ROIs Vector and Rotation Matrix
%
xvec = [1 0 0];         % X-vector for defining 3D cylinder coordinate systems
vec1 = repmat([0 0 1],2,1);  % Cylinders axes is inferior-superior in MRI coordinate system
r1 = eye(3);            % No rotation (MRI coordinate system)
%
% ROI Cylinders
%
[xc,yc,zc] = cylinder(rrad*ones(4,1),72);
zc = 30*zc-15;          % Scale and translate Z coordinates
%
[nrow,ncol] = size(xc);
xyzc = [xc(:) yc(:) zc(:)];
%
% Get Subject Directory
%
% rdirs1 = dir('0?');                    % L: drive directory
% rdirs2 = dir('1?');                    % L: drive directory
% rdirs = [rdirs1; rdirs2];              % L: drive directory
%
rdirs = dir('MRIR*');   % Mack's directory
%
rdirs = {rdirs([rdirs.isdir]').name}';
%
% Select Subject Directory for Creating the Figure
%
idx = listdlg('ListString',rdirs,'Name','Subjects','PromptString', ...
              {'Select subject directory for'; ...
              '       creating the figure.'},'SelectionMode','single');
%
if isempty(idx)
  return;
end
%
rdir1 = rdirs{idx};
%
% Get Visit Subdirectory
%
% rsdirs = {'Visit 1'; 'Visit 2'};       % L: drive directory
rsdirs = {'Visit1'; 'Visit2'};         % Mack's directory
%
idx = menu('Select a visit.','Visit 1','Visit 2');
%
rdir2 = rsdirs{idx};
%
% T1rho Segmentations
%
rdir = fullfile(rdir1,rdir2,'RHO');    % Directory for T1rho segmentations
%
% Get Left or Right Knee
%
leg = ['L'; 'R'];
idx = menu('Select a knee.','Left','Right');
leg = leg(idx);
%
% Get Loaded or Unloaded Knee
%
ld = ['LD'; 'UL'];
idx = menu('Select a knee load.','Loaded','Unloaded');
ld = ld(idx,:);
%
% Read ROIs
%
brois = rd_rois3(rdir,leg,ld,itroch,4);
%
% Get Femur Slices
%
rslc = brois(1).rois(1).slice;         % Femur cartilage
rslb = brois(1).rois(2).slice;         % Femur bone
rslf = intersect(rslc,rslb);           % Ensure unique slices in sorted order
%
% Get Tibia Slices
%
rslc = brois(2).rois(1).slice;         % Tibia cartilage
rslb = brois(2).rois(2).slice;         % Tibia bone
rslt = intersect(rslc,rslb);           % Ensure unique slices in sorted order
rsl = intersect(rslt,rslf);            % Ensure femur and tibia slices
nrsl = size(rsl,1);
%
icmprt = ones(nrsl,1);  % Tibal compartments (1 - lateral, 2 - medial)
%
% Get Cartilage ROI Data for the Femur and Tibia
%
lc = 1;                 % Index to cartilage
sls = cell(2,1);        % Slice numbers for the 1 - femur and 2 - tibia
nsls = zeros(2,1);      % Number of 1 - femur and 2 - tibia slices
%
for lb = 1:2            % Loop through bones - femur = 1 and tibia = 2
%
   sls{lb} = brois(lb).rois(lc).slice; % Slice numbers
   nsls(lb) = length(sls);             % Number of slices
%
   if lb==1
     f3 = cell(nsls(lb),1);            % Femur cartilage coordinates
   else
     t3 = cell(nsls(lb),1);            % Tibia cartilage coordinates
   end
%
   if lb==1&&itroch
     nc = 3;
   else
     nc = 2;
   end
%
   for n = 1:nc         % Loop through compartments - lateral = 1, medial = 2 and trochlea = 3
      [~,~,id] = intersect(brois(lb).rois(lc).roi(n).imageno',sls{lb});
%
      if lb==1                
        f3(id) = brois(lb).rois(lc).roi(n).data3';    % Femur
      else
        t3(id) = brois(lb).rois(lc).roi(n).data3';    % Tibia
        [~,ic] = intersect(rsl,sls{lb}(id));
        icmprt(ic) = n;
      end
   end                  % End of n loop - lateral/medial/trochlear loop
end                     % End of lb loop - femur/tibia loop
%
% Get Indexes to Each Compartment (Lateral/Medial) on Each Bone (Femur/Tibia)
%
ilat = icmprt==1;       % Index to lateral slices
rslc = cell(2,1);       % Slices in 1 - lateral and 2 - medial compartments
rslc{1} = rsl(ilat);    % Lateral slice numbers
rslc{2} = rsl(~ilat);   % Medial slice numbers
%
idxbc = cell(2,2);      % ROI slices (Row: 1 - femur/2 - tibia/Columns: 1 - lateral/2 - medial)
idxplt = cell(2,1);     % Not ROI slices (1 - femur/2 - tibia)
%
for lb = 1:2            % Loop through bones - femur = 1 and tibia = 2
%
   for lc = 1:2         % Loop through compartments - lateral = 1 and medial = 2
%   
      [~,idxbc{lb,lc}] = intersect(sls{lb},rslc{lc},'stable');
%
   end
%
   [~,idxplt{lb}] = setdiff(sls{lb},cell2mat(rslc),'stable');
%
end
%
% Plot Whole Compartment Regions of Interest (ROIs)
%
figure;                 % Knee orientation 
amat = eye(3);
amat(2,2) = -1;         % MRI coordinates to quasi bone coordinates
ax_arrow2(50,1,amat);
view(vw);
axis equal;
axis off;
title('A) Knee Orientation','FontSize',24,'FontWeight','Bold');
%
figure;                 % Combined femur and tibia
for lb = 1:2
   if lb==1
     dat = f3;
   else
     dat = t3;
   end
   plt_datsl(dat(idxplt{lb}),'k',1);
   for lc = 1:2
      plt_datsl(dat(idxbc{lb,lc}),'b',2);
   end
end
view(vw);
axis equal;
axis off;
title('B) Combined Femur and Tibia','FontSize',24,'FontWeight','Bold');
%
ttxt{2} = 'D) Tibia';
ttxt{1} = 'C) Femur';
%
for lb = 1:2
   figure;              % Femur and tibia
   if lb==1
     dat1 = f3;
     dat2 = t3;
   else
     dat1 = t3;
     dat2 = f3;
   end
   plt_datsl(dat2,'k',1);
   plt_datsl(dat1(idxplt{lb}),'k',1);
   for lc = 1:2
      plt_datsl(dat1(idxbc{lb,lc}),'b',2);
   end
   view(vw);
   axis equal;
   axis off;
   title(ttxt{lb},'FontSize',24,'FontWeight','Bold');
end
%
ttxt{4} = 'H) Medial Tibia';
ttxt{3} = 'G) Lateral Tibia';
ttxt{2} = 'F) Medial Femur';
ttxt{1} = 'E) Lateral Femur';
%
for lb = 1:2
%
   if lb==1
     dat1 = f3;
     dat2 = t3;
   else
     dat1 = t3;
     dat2 = f3;
   end
%
   for lc = 1:2
      idt = 2*lb+lc-2;
      figure;                 % Femur and tibia
      plt_datsl(dat2,'k',1);
%      
      plt_datsl(dat1(idxplt{lb}),'k',1);
      plt_datsl(dat1(idxbc{lb,3-lc}),'k',1);
      plt_datsl(dat1(idxbc{lb,lc}),'b',2);
%
      view(vw);
      axis equal;
      axis off;
      title(ttxt{idt},'FontSize',24,'FontWeight','Bold');
   end
%
end
%
% Find Cylindrical Regions of Interest (ROIs)
%
cid = cell(nrsl,1);     % Index to tibia contact points
cdatt = cell(nrsl,1);   % Tibia contact points
cdatf = cell(nrsl,1);   % Femur contact points
%
% Loop through Tibia Slices
%
for k = 1:nrsl          % Loop through slices
%
   slk = rsl(k);        % Slice number
   idf = sls{1}==slk;   % Index to femur slice
   idt = sls{2}==slk;   % Index to tibia slice
%
% Get Contact Points
%
   [cid{k},cdatf{k}] = pwl_contct(f3{idf},t3{idt},thresh);
   if ~isempty(cid{k})
     cdatt{k} = t3{idt}(cid{k},:);
   end
%
end
%
% Get Centers of ROIs from the Tibial Contact Points
%
xyzls = cell2mat(cdatt(ilat));         % Lateral tibia XYZ coordinates in mm
xyzms = cell2mat(cdatt(~ilat));        % Medial tibia XYZ coordinates in mm
%
xyzcc(2,:) = mean(xyzms);              % Center of medial contact ROI
xyzcc(1,:) = mean(xyzls);              % Center of lateral contact ROI
%
% Transform Cartilage Segmentations to ROI Cylinders Axes
%
f3t = cell(2,1);        % Loaded ROIs femur bone transformed (1 - lateral, 2 - medial)
t3t = cell(2,1);        % Loaded ROIs tibia bone transformed (1 - lateral, 2 - medial)
%
for lc = 1:2            % Loop through compartments
   idlat = icmprt==lc;
   [~,idlatf] = intersect(sls{1},rsl(idlat),'stable');
   [~,idlatt] = intersect(sls{2},rsl(idlat),'stable');
   [f3t{lc},t3t{lc}] = coord_tf(xyzcc(lc,:),r1,f3(idlatf),t3(idlatt));
end
%
% Index to Compartment Slices
%
islc = zeros(nrsl,1);
%
ncl = sum(ilat,1);      % Number of lateral slices
islc(ilat) = (1:ncl)';  % Lateral compartment
%
ncm = sum(~ilat,1);     % Number of medial slices
islc(~ilat) = (1:ncm)'; % Medial compartment
%
ipsl = false(ncl,1);    % Index into rslc{1} to contact ROI slices
ipsm = false(ncm,1);    % Index into rslc{2} to contact ROI slices
%
ips = {ipsl; ipsm};     % Lateral and medial combined in cell array
%
l0 = zeros(ncl,1);
m0 = zeros(ncm,1);
nps = {l0 m0; l0 m0};   % Number of segments outside the cylinder
%
% Initialize Arrays
%
xyzpi = {cell(ncl,1) cell(ncm,1); cell(ncl,1) cell(ncm,1)};
xyzpo = {cell(ncl,2) cell(ncm,2); cell(ncl,2) cell(ncm,2)};
%
% Loop through the Slices to Define Analysis Slices
%
for k = 1:nrsl          % Loop through slices
%
% Find Lateral and Medial Loaded ROIs on Slice
%
   slk = rsl(k);        % Slice number
   l = icmprt(k);       % Compartments
   ks = islc(k);        % Index to compartment slice
%
% Get Coordinates for Lines
%
   fb = f3t{l}{ks};     % Slice segmentations in ROI coordinates
   tb = t3t{l}{ks};     % Slice segmentations in ROI coordinates
%
% Get ROIs in Rotated Coordinates Based on Intersections with
% Cylindrical ROIs
%
   [xyzfi,xyzfo] = cyl_pwl1(fb,rrad,15);
   [xyzti,xyzto] = cyl_pwl1(tb,rrad,15);
%
% Plot Slice
%
   if iplt
     figure;
     plot3(fb(:,1),fb(:,2),fb(:,3),'g.-');
     hold on;
     plot3(tb(:,1),tb(:,2),tb(:,3),'b.-');
   end
%
% Intersections with the Loaded Cylindrical ROIs?
%
   ifemur = ~isempty(xyzfi);
   itibia = ~isempty(xyzti);
%
   if ifemur
%
     ips{l}(ks) = true; % Contact ROI on slice
     nps{1,l}(ks) = size(xyzfo,1);     % Number of segments outside cylinder
%
     if iplt
       plot3(xyzfi(:,1),xyzfi(:,2),xyzfi(:,3),'rs');
     end
%
     [xyzif,xyzof] = tf_coord(r1,xyzcc(l,:),{xyzfi},xyzfo);     % Transform Back to MRI Coordinates
%
     xyzpi{1,l}(ks) = xyzif;           % Put Segments into Cell Arrays
     xyzpo{1,l}(ks,1:nps{1,l}(ks)) = xyzof';     % Put Segments into Cell Arrays
%
   end
%
   if itibia
%
     ips{l}(ks) = true; % Contact ROI on slice
     nps{2,l}(ks) = size(xyzto,1);     % Number of segments outside cylinder
%
     if iplt
       plot3(xyzti(:,1),xyzti(:,2),xyzti(:,3),'rs');
     end
%
     [xyzit,xyzot] = tf_coord(r1,xyzcc(l,:),{xyzti},xyzto);     % Transform Back to MRI Coordinates
%
     xyzpi{2,l}(ks) = xyzit;           % Put Segments into Cell Arrays
     xyzpo{2,l}(ks,1:nps{2,l}(ks)) = xyzot';     % Put Segments into Cell Arrays
%
   end
%
   if (ifemur||itibia)&&iplt
     pause;
   end
%
end
%
% Get Slices with Contact Regions of Interest (ROIs)
%
lidx = cellfun(@logical,nps,'UniformOutput',false);   % Logical index to slices
pslc = [rslc'; rslc'];  % Compartment slice numbers
%
idxcplt = cell(2,2);    % Not ROI slices (Row: 1 - femur/2 - tibia/Columns: 1 - lateral/2 - medial)
%
% Loop through Contact Regions of Interest (ROIs)
%
ttxt{4} = 'L) Medial Tibia Contact';
ttxt{3} = 'K) Lateral Tibia Contact';
ttxt{2} = 'J) Medial Femur Contact';
ttxt{1} = 'I) Lateral Femur Contact';
%
for lb = 1:2            % lb = 1 - femur/lb = 2 - tibia
%
   if lb==1
     dat1 = f3;
     dat2 = t3;
   else
     dat1 = t3;
     dat2 = f3;
   end
%
   for lc = 1:2         % lc = 1 - lateral/lc = 2 - medial
      idt = 2*lb+lc-2;
%
      nps{lb,lc} = nps{lb,lc}(lidx{lb,lc});
      pslc{lb,lc} = pslc{lb,lc}(lidx{lb,lc});
      xyzpi{lb,lc} = xyzpi{lb,lc}(lidx{lb,lc});
      xyzpo{lb,lc} = [xyzpo{lb,lc}(lidx{lb,lc},1) ...
                      xyzpo{lb,lc}(lidx{lb,lc},2)];

      [~,idxcplt{lb,lc}] = setdiff(sls{lb},pslc{lb,lc},'stable');
%
      figure;                 % Both femur and tibia
      plt_datsl(dat2,'k',1);
%
      plt_datsl(dat1(idxcplt{lb,lc}),'k',1);
      idxsl = nps{lb,lc}>=1;
      plt_datsl(xyzpo{lb,lc}(idxsl,1),'k',1);
      idxsl = nps{lb,lc}==2;
      plt_datsl(xyzpo{lb,lc}(idxsl,2),'k',1);
      plt_datsl(xyzpi{lb,lc},'b',2);
%
      view(120,37.5);
      axis equal;
      axis off;
      title(ttxt{idt},'FontSize',24,'FontWeight','Bold');
   end
%
end
%
return