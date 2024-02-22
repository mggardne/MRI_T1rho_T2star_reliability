%#######################################################################
%
%    * SEGmentation to Two (2) Regions of Interest (ROIs) Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     CSV files to create masks for cylindrical regions of interest in
%     the lateral and medial tibial compartments and corresponding
%     femoral ROIs.  In addition, masks are created for unloaded
%     cylindrical regions of interest in the lateral and medial
%     posterior femoral condyles. The masks are saved in MAT files with
%     the series number and ending in "_2rois.mat."
%
%     The center of the tibial cylindrical ROIs is determined using the
%     center of the contact area on the tibia and femur cartilage
%     surface.
%
%     The center of the unloaded femoral cylindrical ROIs is determined
%     as the center of the unloaded posterior condyles.  See CSV file,
%     UnloadedFemPlugs.csv.
%
%     NOTES:  1.  The registered MRI MAT files must be in subject
%             directories starting with "MRIR" and either "Visit1" or
%             "Visit2" visit subdirectories.
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  See rd_dicom.m.
%
%             3.  The M-files coord_tf.m, cr_mask2.m, cr_mask2f.m,
%             cyl_lin3.m, cyl_pl3.m, cyl_pwl.m, decomp.m, fix_gap.m,
%             in_tri2d.m, lsect2.m, lsect2a.m, lsect3.m, lsect4.m,
%             lsect5.m, midline.m, mk2_tri_2d.m, mk2_tri_2df.m,
%             pts2lin.m, mk2_tri_2d.m, pts2lin.m, pwl_contct.m,
%             rd_rois3.m, rd_roi6.m, tf_coord.m, trnsf2pixel.m, and
%             vec2rot3d.m must be in the current directory or path.
%
%             4.  CSV file, UnloadedFemPlugs.csv, must be in the main
%             MRI reliability study folder, "MRI_Reliability_Study\"
%             (two directories above the visit subdirectories).
%
%     20-Sep-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Setup Color Maps, Line Types and Color, and Legend Strings
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
rmap = [[0.85 0 0]; [0 0 0.85]; gmap]; % Red/blue color map for ROIs
%
lt = ['g.-'; 'b.-'; 'c.-'; 'r.-'; 'y.-'; 'm.-']; % Line color and type
legds = ['FEM_CART'; 'FEM_BONE'; 'TIB_CART'; 'TIB_BONE'];
%
% Parameter for Dividing Cartilage in Half and Radius of Analysis Region
%
dist = 7.5;             % Maximum distance to midline in pixels
rrad = 6.0;             % Radius of analysis region (width of 3 or 4 slices)
rrad2 = rrad*rrad;
itroch = true;          % Read trochlea?
% itroch = false;         % Read trochlea?
%
% Axis Vector and Contact ROIs Vector and Rotation Matrix
%
xvec = [1 0 0];         % X-vector for defining 3D cylinder coordinate systems
vec1 = repmat([0 0 1],2,1);  % Cylinders axes is inferior-superior in MRI coordinate system
% r1 = eye(3);            % No rotation (MRI coordinate system)
%
% ROI Cylinders
%
[xc,yc,zc] = cylinder(rrad*ones(4,1),72);
zc = 30*zc-15;          % Scale and translate Z coordinates
%
[nrow,ncol] = size(xc);
xyzc = [xc(:) yc(:) zc(:)];
%
% Read Unloaded Femoral Condyles ROIs Centers and Vectors from CSV File:
% UnloadedFemPlus.csv
%
dat = readmatrix(fullfile('..','..','UnloadedFemPlugs.csv'), ...
                 'NumHeaderLines',1);
fidx = dat(:,1:6);      % Indices
fpts = dat(:,8:10);     % Center points coordinates
fvec = dat(:,11:13);    % ROIs normals (vectors)
fvec = fvec./sqrt(sum(fvec.*fvec,2));  % Normalize vectors
%
fsubj = fidx(:,1);      % Subject number
fvisit = fidx(:,2)+1;   % Visit number
fres = fidx(:,3);       % T1rho = 0 and T2* = 1
fleg = fidx(:,4);       % Leg
fload = fidx(:,5);      % Load
fcmprt = fidx(:,6)+1;   % Compartment
%
clear dat fidx;
%
% Compartment Labels
%
cmprt = {'Lateral'; 'Medial'};
%
% Directory with MAT Files as a String
%
dirstr = split(pwd,filesep);
vstr = dirstr{end};
vstr = [vstr(1:end-1) ' ' vstr(end)];
dirstr = dirstr{end-1};
dirstr = split(dirstr,' ');
dirstr = [dirstr{end} ' - ' vstr];
%
% Get Subject and Visit Indices
%
subj = eval(dirstr(1:2));
fidx = find(fsubj==subj);
%
vnum = eval(dirstr(end));
fv = find(fvisit==vnum);
%
fidx = intersect(fidx,fv);
%
% Skip T1rho?
%
% iskip = true;
iskip = false;
%
if ~iskip
%
% T1rho
%
rdir = 'RHO';           % Directory for T1rho segmentations
id5m = 1;               % Use first spin lock time for plots
thresh = 1.953125;      % Four (4) pixels = 4*0.48828125 = 1.953125
% threshs = num2str(thresh,'%.3f');      % Threshold as a string
fr = find(fres==0);
fidr = intersect(fidx,fr);             % Subject, visit and T1rho analysis
%
% Get T1rho MAT Files in Directory
%
d = dir('T1rho_S*.mat');
mnams = {d.name}';
idr = contains(mnams,'roi','IgnoreCase',true);
mnams = mnams(~idr);
id3 = contains(mnams,'_3_');
mnams = mnams(~id3);
nmat = size(mnams,1);
%
if nmat~=4
  error(' *** ERROR in seg_2rois:  Not four (4) MAT files!')
end
%
% Loop through T1rho MAT Files
%
% for m = 1:nmat
for m = 4
%
   mnam = mnams{m};
   load(mnam);
   fs = ['S' snt];      % Series number prefaced with a 'S'
%
% Parse Series Text for Leg and Load
%
   if strcmpi(st(1),'L')
     leg = 'L';
     fl = find(fleg==0);
   else
     leg = 'R';
     fl = find(fleg==1);
   end
%
   fids = intersect(fidr,fl);          % Include index to leg
%
   if contains(st,'Load','IgnoreCase',true)
     ld = 'LD';
     fl = find(fload==1);              % Loaded
   else
     ld = 'UL';
     fl = find(fload==0);              % Unloaded
   end
%
   fids = intersect(fids,fl);          % Include index to load
%
% Read ROIs
%
   brois = rd_rois3(rdir,leg,ld,itroch,4);
%
% Get Femur Slices
%
   rslc = brois(1).rois(1).slice;      % Femur cartilage
   rslb = brois(1).rois(2).slice;      % Femur bone
   rslf = intersect(rslc,rslb);        % Ensure unique slices in sorted order
%
% Get Tibia Slices
%
   rslc = brois(2).rois(1).slice;      % Tibia cartilage
   rslb = brois(2).rois(2).slice;      % Tibia bone
   rsl = intersect(rslc,rslb);         % Ensure unique slices in sorted order
   rsl = intersect(rsl,rslf);          % Ensure femur and tibia slices
   nrsl = size(rsl,1);
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_2ROIs1.ps'];          % ROI lines print file name
   pnam2 = [fs '_2ROIs2.ps'];          % ROI areas print file name
   pnam3 = [fs '_2ROIs3.ps'];          % ROI regions cylinders print file name
   pnam4 = [fs '_2ROIs4.ps'];          % ROI regions masks print file name
%
   f = cell(2,nrsl);    % Femur coordinates (1 - cartilage, 2 - bone)
   t = cell(2,nrsl);    % Tibia coordinates (1 - cartilage, 2 - bone)
%
   f3 = cell(2,nrsl);   % Femur coordinates (1 - cartilage, 2 - bone)
   t3 = cell(2,nrsl);   % Tibia coordinates (1 - cartilage, 2 - bone)
%
   ibone = false(nrsl,2);              % Bone exists? (1 - femur, 2 - tibia)
   icmprt = ones(nrsl,1);              % Tibal compartments (1 - lateral, 2 - medial)
   npxt = prod(iszs);                  % Total number of pixels in a slice
   maskf = false(npxt,2,nrsl);         % Mask for femoral cartilage
   maskt = false(npxt,2,nrsl);         % Mask for tibial cartilage
%
   cid = cell(nrsl,1);                 % Index to tibia contact points
   cdatt = cell(nrsl,1);               % Tibia contact points
   cdatf = cell(nrsl,1);               % Femur contact points
%
   xyz2 = zeros(2,3);   % Unloaded posterior femoral center points
   vec2 = zeros(2,3);   % Unloaded posterior femoral normal vectors
%
% Loop through Tibia Slices
%
%    for k = 1:nrsl       % Loop through slices
   for k = 6:12         % Loop through slices
%
% Get T1 Slice Image
%
      slk = rsl(k);        % Slice number
      img = squeeze(v(:,:,slk,id5m));
%
      figure;
      orient landscape;
      imagesc(img);
      colormap gray;
      axis image;
      axis off;
      title({dirstr; [fs ' Slice ' int2str(slk)]},'FontSize',16, ...
            'FontWeight','bold');
      hold on;
%
      lh = gobjects(4,1);              % Line graphic handles
      idl = false(4,1);
%
% Get ROI Data for this Slice and Plot ROIs
%
      for lb = 1:2      % Loop through bones - femur = 1 and tibia = 2
         for lc = 1:2   % Loop through surfaces - cartilage = 1 and bone = 2
            idxl = 2*lb+lc-2;
            idxr = brois(lb).rois(lc).slice==slk;
            if any(idxr)
              if lb==1&&itroch
                nc = 3;
              else
                nc = 2;
              end
              for n = 1:nc   % Loop through compartments - lateral = 1, medial = 2 and trochlea = 3
                 idxs = brois(lb).rois(lc).roi(n).imageno==slk;
                 if any(idxs)
                   dat = cell2mat( ...
                                 brois(lb).rois(lc).roi(n).data(idxs)');
                   dat3 = cell2mat( ...
                                brois(lb).rois(lc).roi(n).data3(idxs)');
                   [dat,idx] = fix_gap(dat);
                   dat3 = dat3(idx,:);
                   if lb==1
                     f{lc,k} = dat;    % Femur
                     f3{lc,k} = dat3;  % Femur
                   else
                     t{lc,k} = dat;    % Tibia
                     t3{lc,k} = dat3;  % Tibia
                     icmprt(k) = n;
                   end
                   lh(idxl) = plot(dat(:,1),dat(:,2),lt(idxl,:));
                   idl(idxl) = true;
                 end
              end       % End of n loop - lateral/medial loop
            end         % End of if - any rois for this slice
         end            % End of lc loop - cartilage/bone loop
         if lb==1
           ibone(k,lb) = all(~cellfun('isempty',f(:,k)));
         else
           ibone(k,lb) = all(~cellfun('isempty',t(:,k)));
         end
      end               % End of lb loop - femur/tibia loop
%
% Get Transform Parameters
% Scale, s, Rotation Matrix, r, and Translation Vector, tv.
%
      if k==1
        [~,~,r,tv] = trnsf2pixel(f3{2,k},f{2,k},f3{2,k});
        s = 1./pspcs(1);               % Scale
        r3 = -r(:,[1 3 2]);
      end
%
% Get Contact Points
%
      [cid{k},cdatf{k}] = pwl_contct(f3{1,k},t3{1,k},thresh);
      if ~isempty(cid{k})
        cdatt{k} = t3{1,k}(cid{k},:);
      end
%
% Add Legends and Print Slice Plots
%
      legend(lh(idl),legds(idl,:),'Interpreter','none');
      k
      pause
%       if k==1
%         print('-dpsc2','-r600','-fillpage',pnam1);
%       else
%         print('-dpsc2','-r600','-fillpage','-append',pnam1);
%       end
%
% Create Logical Masks for the Cartilage on this Slice
%
      if ibone(k,1)
        [maskf(:,1,k),maskf(:,2,k)] = cr_mask2f(f(:,k),iszs,dist,pspcs);
      end
%
      if ibone(k,2)
        [maskt(:,1,k),maskt(:,2,k)] = cr_mask2(t(:,k),iszs,dist,pspcs);
      end
%
% Plot ROIs
%
      mask1 = maskf(:,1,k)|maskt(:,1,k);    % Superficial cartilage mask
      mask2 = maskf(:,2,k)|maskt(:,2,k);    % Deep cartilage mask
      mask = mask1|mask2;              % All cartilage mask
      cmx = max(img(:));
      img1 = img-cmx-1;
      idmsk = ~mask;
      img2 = img;
      img2(idmsk) = img(idmsk)-max(img(:))-1;
      dcmx = 16*cmx/128;
      img1(mask1) = dcmx;  % Blue - Superficial
      img1(mask2) = cmx-dcmx;   % Red - Deep
%
      figure;
      orient landscape;
      imagesc(img1);
      colormap(cmap);
      caxis([-cmx cmx]);
      axis image;
      axis off;
      title({dirstr; [fs ' Slice ' int2str(slk)]},'FontSize',16, ...
            'FontWeight','bold');
%
%       if k==1
%         print('-dpsc2','-r600','-fillpage',pnam2);
%       else
%         print('-dpsc2','-r600','-fillpage','-append',pnam2);
%       end
%
   end                  % End of k loop - tibia slices
%
return
%
% Close Slice Plots
%
%    close all;
%
% Get Centers of ROIs from the Tibial Contact Points
%
   ilat = icmprt==1;    % Lateral compartment
   xyzls = cell2mat(cdatt(ilat));      % Lateral tibia XYZ coordinates in mm
   xyzms = cell2mat(cdatt(~ilat));     % Medial tibia XYZ coordinates in mm
%
%    idv = ~cellfun('isempty',cdatt);
%    nt = cellfun('size',cdatt,1);       % Number of contact points on each slice
%
   xyz1(2,:) = mean(xyzms);            % Center of medial contact ROI
   xyz1(1,:) = mean(xyzls);            % Center of lateral contact ROI
%
% Get Centers and Axial Vectors of Posterior Femoral ROIs
% Row 1 -> Lateral, and Row 2 -> Medial
%
   idc2 = fcmprt(fids);
   xyz2(idc2,:) = fpts(fids,:);        % Centers of cylinders
   vec2(idc2,:) = fvec(fids,:);        % Cylinder axes
%
% Transform Cartilage and Femoral Segmentations to ROI Cylinders Axes
%
   ft1 = cell(2,1);     % Loaded ROIs femur-tibia bone transformed (1 - lateral, 2 - medial)
   r2 = zeros(3,3,2);   % Rotation matrix to cylinders axes
   ft2 = cell(2,1);     % Unloaded ROIs femur transformed (1 - lateral, 2 - medial)
%
   for l = 1:2          % Loop through compartments
      ilat = icmprt==l;
%       ft1{l} = coord_tf(xyz1(l,:),r1,[f3(2,ilat); t3(2,ilat)]);
      ft1{l} = coord_tf(xyz1(l,:),r3,[f3(2,ilat); t3(2,ilat)]);
      r2(:,:,l) = vec2rot3d(vec2(l,:),xvec);
      ft2{l} = coord_tf(xyz2(l,:),r2(:,:,l),f3(:,ilat));
   end
%
% Index to Compartment Slices
%
   islc = zeros(nrsl,1);
%
   idcl = find(icmprt==1);             % Lateral compartment
   ncl = size(idcl,1);
   islc(idcl) = (1:ncl)';
%
   idcl = find(icmprt==2);             % Medial compartment
   ncl = size(idcl,1);
   islc(idcl) = (1:ncl)';
%
%  Setup Figures
%
   hf3 = gobjects(5,1);  % Figure handles for rotated cylindrical ROIs
   hf3(5) = figure;      % Figure handle for cylindrical ROIs
   orient landscape;
%
% Plot Cylinders in MRI Coordinates
%
%    xyzct = tf_coord(r1,xyz1(1,:),{xyzc});        % Loaded lateral ROI
   xyzct = tf_coord(r3',xyz1(1,:),{xyzc});       % Loaded lateral ROI
   xyzct = xyzct{1};
   xct = reshape(xyzct(:,1),nrow,ncol);
   yct = reshape(xyzct(:,2),nrow,ncol);
   zct = reshape(xyzct(:,3),nrow,ncol);
%
   plot3(xct,yct,zct,'y-','LineWidth',1);
   hold on;
   plot3(xct',yct',zct','y-','LineWidth',1);
   axis equal;
%
%    xyzct = tf_coord(r1,xyz1(2,:),{xyzc});        % Loaded medial ROI
   xyzct = tf_coord(r3',xyz1(2,:),{xyzc});       % Loaded medial ROI
   xyzct = xyzct{1};
   xct = reshape(xyzct(:,1),nrow,ncol);
   yct = reshape(xyzct(:,2),nrow,ncol);
   zct = reshape(xyzct(:,3),nrow,ncol);
%
   plot3(xct,yct,zct,'y-','LineWidth',1);
   plot3(xct',yct',zct','y-','LineWidth',1);
%
   xyzct = tf_coord(r2(:,:,1)',xyz2(1,:),{xyzc});     % Unloaded lateral ROI
   xyzct = xyzct{1};
   xct = reshape(xyzct(:,1),nrow,ncol);
   yct = reshape(xyzct(:,2),nrow,ncol);
   zct = reshape(xyzct(:,3),nrow,ncol);
%
   plot3(xct,yct,zct,'y-','LineWidth',1);
   plot3(xct',yct',zct','y-','LineWidth',1);
%
   xyzct = tf_coord(r2(:,:,2)',xyz2(2,:),{xyzc});     % Unloaded medial ROI
   xyzct = xyzct{1};
   xct = reshape(xyzct(:,1),nrow,ncol);
   yct = reshape(xyzct(:,2),nrow,ncol);
   zct = reshape(xyzct(:,3),nrow,ncol);
%
   plot3(xct,yct,zct,'y-','LineWidth',1);
   plot3(xct',yct',zct','y-','LineWidth',1);
%
% Plot Cylinders in Rotated Coordinates Aligned with the Axis of the
% Cylinders
%
   hf3(1) = figure;     % Lateral loaded ROI
   orient landscape;
%
   plot3(xc,yc,zc,'y-');
   hold on;
   plot3(xc',yc',zc','y-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Lateral Loaded ROI','FontSize',14,'FontWeight','bold');
%
   hf3(2) = figure;     % Medial loaded ROI
   orient landscape;
%
   plot3(xc,yc,zc,'y-');
   hold on;
   plot3(xc',yc',zc','y-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Medial Loaded ROI','FontSize',14,'FontWeight','bold');
%
   hf3(3) = figure;     % Lateral unloaded ROI
   orient landscape;
%
   plot3(xc,yc,zc,'y-');
   hold on;
   plot3(xc',yc',zc','y-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Lateral Unloaded ROI','FontSize',14,'FontWeight','bold');
%
   hf3(4) = figure;     % Medial unloaded ROI
   orient landscape;
%
   plot3(xc,yc,zc,'y-');
   hold on;
   plot3(xc',yc',zc','y-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Medial Unloaded ROI','FontSize',14,'FontWeight','bold');
%
% Loop through the Slices to Define Analysis Regions
%
   irsl1 = false(nrsl,2);    % Loaded contact ROIs on slice (column 1 - lateral, column 2 - medial)
   irsl2 = false(nrsl,2);    % Unloaded posterior ROIs on slice
%
   maskr1 = false(npxt,2,nrsl);        % Mask for loaded ROI
   maskr2 = false(npxt,2,nrsl);        % Mask for unloaded ROI
%
   for k = 1:nrsl       % Loop through slices
%
% Plot Unrotated Segmentations
%
      figure(hf3(5));
%
      plot3(f3{1,k}(:,1),f3{1,k}(:,2),f3{1,k}(:,3),lt(1,:), ...
            'LineWidth',1);
      plot3(f3{2,k}(:,1),f3{2,k}(:,2),f3{2,k}(:,3),lt(2,:), ...
            'LineWidth',1);
      plot3(t3{1,k}(:,1),t3{1,k}(:,2),t3{1,k}(:,3),lt(3,:), ...
            'LineWidth',1);
      plot3(t3{2,k}(:,1),t3{2,k}(:,2),t3{2,k}(:,3),lt(4,:), ...
            'LineWidth',1);
%
% Find Lateral and Medial Loaded and Unloaded ROIs on Slice
%
      l = icmprt(k);    % Compartments
      ks = islc(k);
%
% Get Coordinates for Both Lines
%
      fb1 = ft1{l}{1,ks};    % Slice segmentations in loaded ROI coordinates
      tb1 = ft1{l}{2,ks};
%
      fb2 = ft2{l}{1,ks};    % Slice segmentations in unloaded ROI coordinates
      tb2 = ft2{l}{2,ks};
%
% Get Polygon ROIs in Rotated Coordinates Based on Intersections with
% Cylindrical ROIs
%
      [xyzp1,xyzif1,xyzit1] = cyl_pwl(fb1,tb1,rrad,15,24);
      [xyzp2,xyzif2,xyzit2] = cyl_pwl(fb2,tb2,rrad,15,12);
%
% Intersections with the Loaded Cylindrical ROIs?
%
      if ~isempty(xyzp1)
%
        irsl1(k,l) = true;   % Intersection true for this slice and compartment
%
% Plot Lines and Polygon
%
        figure(hf3(l));
%
        plot3(fb1(:,1),fb1(:,2),fb1(:,3),'b.-');
        plot3(tb1(:,1),tb1(:,2),tb1(:,3),'g.-');
%
        plot3(fb1(1,1),fb1(1,2),fb1(1,3),'b^');            % Start
        plot3(fb1(end,1),fb1(end,2),fb1(end,3),'bs');      % End
%
        plot3(tb1(1,1),tb1(1,2),tb1(1,3),'g^');            % Start
        plot3(tb1(end,1),tb1(end,2),tb1(end,3),'gs');      % End
%
        if ~isempty(xyzif1)
          plot3(xyzif1(:,1),xyzif1(:,2),xyzif1(:,3),'rs');
        end
        if ~isempty(xyzit1)
          plot3(xyzit1(:,1),xyzit1(:,2),xyzit1(:,3),'rs');
        end
%
        hp = patch(xyzp1(:,1),xyzp1(:,2),xyzp1(:,3),'k');
        set(hp,'FaceAlpha',0.4);
        set(hp,'EdgeColor','none');
%
% Transform Back to MRI Coordinates
%
%         xyzpt1 = tf_coord(r1,xyz1(l,:),{xyzp1});
        xyzpt1 = tf_coord(r3',xyz1(l,:),{xyzp1});
        xyzpt1 = xyzpt1{1};
%
        figure(hf3(5));
        hp = patch(xyzpt1(:,1),xyzpt1(:,2),xyzpt1(:,3),'k');
        set(hp,'FaceAlpha',0.4);
        set(hp,'EdgeColor','none');
%
% Transform to Pixel Coordinates
%
        xyzpp1 = s*xyzpt1*r+repmat(tv,size(xyzpt1,1),1);
%
% Create Mask from Polygon
%
        msk = poly2mask(xyzpp1(:,1),xyzpp1(:,2),iszs(1),iszs(2));
        maskr1(:,l,k) = msk(:);
%
      end
%
% Intersections with the Unloaded Cylindrical ROIs?
%
      if ~isempty(xyzp2)
%
        irsl2(k,l) = true;   % Intersection true for this slice and compartment
%
% Plot Lines and Polygon
%
        figure(hf3(l+2));
%
        plot3(fb2(:,1),fb2(:,2),fb2(:,3),'b.-');
        plot3(tb2(:,1),tb2(:,2),tb2(:,3),'g.-');
%
        plot3(fb2(1,1),fb2(1,2),fb2(1,3),'b^');            % Start
        plot3(fb2(end,1),fb2(end,2),fb2(end,3),'bs');      % End
%
        plot3(tb2(1,1),tb2(1,2),tb2(1,3),'g^');            % Start
        plot3(tb2(end,1),tb2(end,2),tb2(end,3),'gs');      % End
%
        if ~isempty(xyzif2)
          plot3(xyzif2(:,1),xyzif2(:,2),xyzif2(:,3),'rs');
        end
        if ~isempty(xyzit2)
          plot3(xyzit2(:,1),xyzit2(:,2),xyzit2(:,3),'rs');
        end
%
        hp = patch(xyzp2(:,1),xyzp2(:,2),xyzp2(:,3),'k');
        set(hp,'FaceAlpha',0.4);
        set(hp,'EdgeColor','none');
%
% Transform Back to MRI Coordinates
%
        xyzpt2 = tf_coord(r2(:,:,l)',xyz2(l,:),{xyzp2});
        xyzpt2 = xyzpt2{1};
%
        figure(hf3(5));
        hp = patch(xyzpt2(:,1),xyzpt2(:,2),xyzpt2(:,3),'k');
        set(hp,'FaceAlpha',0.4);
        set(hp,'EdgeColor','none');
%
% Transform to Pixel Coordinates
%
        xyzpp2 = s*xyzpt2*r+repmat(tv,size(xyzpt2,1),1);
%
% Create Mask from Polygon
%
        msk = poly2mask(xyzpp2(:,1),xyzpp2(:,2),iszs(1),iszs(2));
        maskr2(:,l,k) = msk(:);
%
      end
%
   end                  % End of slice loop - k
%
% Print Plots
%
   figure(hf3(5));
   title({dirstr; fs; 'T1\rho ROIs'},'FontSize',12,'FontWeight','bold');
%
   for k = 1:5
%
      figure(hf3(k));
%
      if k==1
        print('-dpsc2','-r600','-fillpage',pnam3);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam3);
      end
%
      if k==5
        view(2);
        print('-dpsc2','-r600','-fillpage','-append',pnam3);
      end
%
   end
%
   close all;
%
% Combine Cartilage and ROI Masks to Get Lateral Loaded Mask
%
   rsll1 = rsl(irsl1(:,1));
   nrsll1 = size(rsll1,1);
%
   maskfrl1 = false(npxt,2,nrsll1);    % Mask for loaded region femoral cartilage
   masktrl1 = false(npxt,2,nrsll1);    % Mask for loaded region tibial cartilage
%
   for k = 1:nrsll1
%
      ida = ismember(rsl,rsll1(k));
%
      for l = 1:2       % 1 - Superficial, 2 - Deep
%
         maskfrl1(:,l,k) = maskf(:,l,ida)&maskr1(:,1,ida);
         masktrl1(:,l,k) = maskt(:,l,ida)&maskr1(:,1,ida);
%
      end
%
   end
%
% Combine Cartilage and ROI Masks to Get Medial Loaded Mask
%
   rslm1 = rsl(irsl1(:,2));
   nrslm1 = size(rslm1,1);
%
   maskfrm1 = false(npxt,2,nrslm1);    % Mask for loaded region femoral cartilage
   masktrm1 = false(npxt,2,nrslm1);    % Mask for loaded region tibial cartilage
%
   for k = 1:nrslm1
%
      ida = ismember(rsl,rslm1(k));
%
      for l = 1:2       % 1 - Superficial, 2 - Deep
%
         maskfrm1(:,l,k) = maskf(:,l,ida)&maskr1(:,2,ida);
         masktrm1(:,l,k) = maskt(:,l,ida)&maskr1(:,2,ida);
%
      end
%
   end
%
% Combine Cartilage and ROI Masks to Get Lateral Unloaded Femoral Mask
%
   rsll2 = rsl(irsl2(:,1));
   nrsll2 = size(rsll2,1);
%
   maskfrl2 = false(npxt,2,nrsll2);    % Mask for unloaded region femoral cartilage
%
   for k = 1:nrsll2
%
      ida = ismember(rsl,rsll2(k));
%
      for l = 1:2       % 1 - Superficial, 2 - Deep
%
         maskfrl2(:,l,k) = maskf(:,l,ida)&maskr2(:,1,ida);
%
      end
%
   end
%
% Combine Cartilage and ROI Masks to Get Medial Unloaded Femoral Mask
%
   rslm2 = rsl(irsl2(:,2));
   nrslm2 = size(rslm2,1);
%
   maskfrm2 = false(npxt,2,nrslm2);    % Mask for unloaded region femoral cartilage
%
   for k = 1:nrslm2
%
      ida = ismember(rsl,rslm2(k));
%
      for l = 1:2       % 1 - Superficial, 2 - Deep
%
         maskfrm2(:,l,k) = maskf(:,l,ida)&maskr2(:,2,ida);
%
      end
%
   end
%
% Plot the ROI Masks
%
   rslp = unique([rsll1; rsll2; rslm1; rslm2]);
   nrslp = size(rslp,1);
%
   for k = 1:nrslp
%
      slk = rslp(k);
      idc = rsl==slk;
      idc = icmprt(idc);
%
      img = squeeze(v(:,:,slk,id5m));
      cmx = max(img(:));
      img1 = img-cmx-1;
%
      figure;
      orient landscape;
%
% Plot Loaded Contact and Unloaded Posterior ROIs
%
      mask1 = false(npxt,1);
      mask2 = false(npxt,1);
%
      idp = rsll1==slk;
      if any(idp)
        mask1 = maskfrl1(:,1,idp)|masktrl1(:,1,idp)|mask1; % Superficial cartilage mask
        mask2 = maskfrl1(:,2,idp)|masktrl1(:,2,idp)|mask2; % Deep cartilage mask
      end
%
      idp = rslm1==slk;
      if any(idp)
        mask1 = maskfrm1(:,1,idp)|masktrm1(:,1,idp)|mask1; % Superficial cartilage mask
        mask2 = maskfrm1(:,2,idp)|masktrm1(:,2,idp)|mask2; % Deep cartilage mask
      end
%
      idp = rsll2==slk;
      if any(idp)
        mask1 = maskfrl2(:,1,idp)|mask1;    % Superficial cartilage mask
        mask2 = maskfrl2(:,2,idp)|mask2;    % Deep cartilage mask
      end
%
      idp = rslm2==slk;
      if any(idp)
        mask1 = maskfrm2(:,1,idp)|mask1;    % Superficial cartilage mask
        mask2 = maskfrm2(:,2,idp)|mask2;    % Deep cartilage mask
      end
%
      dcmx = 16*cmx/128;
      img1(mask1) = dcmx;         % Blue - Superficial
      img1(mask2) = cmx-dcmx;     % Red - Deep
%
      imagesc(img1);
      colormap(cmap);
      caxis([-cmx cmx]);
      axis image;
      axis off;
      title({[fs ' Slice ' int2str(slk)]; [cmprt{idc}, ...
            ' Compartment']; 'T1\rho ROIs'}, ...
            'FontSize',12,'FontWeight','bold');
%
      if k==1
        print('-dpsc2','-r600','-fillpage',pnam4);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam4);
      end
%
   end
%
% Save Masks, ROIS and Slice Information into MAT File
%
   savnam = [mnam(1:end-4) '_2rois.mat'];
   save(savnam,'f','f3','ibone','icmprt','maskf','maskfrl1', ...
               'maskfrm1','maskfrl2','maskfrm2','maskr1','maskr2', ...
               'maskt','masktrl1','masktrm1','brois','rsl', ...
               'rsll1','rslm1','rsll2','rslm2','t','t3');
%
   close all;
%
end                     % End of m loop - MAT file loop
%
end                     % End of if skip RHO
%
% T2star
%
rdir = 'T2S';           % Directory for T2* segmentations
id5 = repmat(3,4,1);    % Default echo time for segmentations
thresh = 1.75;          % Four (4) pixels = 4*0.4375 = 1.75
% threshs = num2str(thresh,'%.3f');      % Threshold as a string
fr = find(fres==1);
fidr = intersect(fidx,fr);             % Subject, visit and T2* analysis
%
% Get T2* MAT Files in Directory
%
d = dir('T2star_S*.mat');
mnams = {d.name}';
idr = contains(mnams,'roi','IgnoreCase',true);
mnams = mnams(~idr);
id3 = contains(mnams,'_3_');
mnams = mnams(~id3);
nmat = size(mnams,1);
%
if nmat~=4
  error(' *** ERROR in seg_2rois:  Not four (4) MAT files!')
end
%
% Loop through T2* MAT Files
%
for m = 1:nmat
%
   mnam = mnams{m};
   load(mnam);
   fs = ['S' snt];      % Series number prefaced with a 'S'
   if length(id5)>=m
     id5m = id5(m);
   else
     id5m = id5(1);
   end
%
% Parse Series Text for Leg and Load
%
   if strcmpi(st(1),'L')
     leg = 'L';
     fl = find(fleg==0);
   else
     leg = 'R';
     fl = find(fleg==1);
   end
%
   fids = intersect(fidr,fl);          % Include index to leg
%
   if contains(st,'Load','IgnoreCase',true)
     ld = 'LD';
     fl = find(fload==1);              % Loaded
   else
     ld = 'UL';
     fl = find(fload==0);              % Unloaded
   end
%
   fids = intersect(fids,fl);          % Include index to load
%
% Read ROIs
%
   brois = rd_rois3(rdir,leg,ld,itroch,5);
%
% Get Femur Slices
%
   rslc = brois(1).rois(1).slice;      % Femur cartilage
   rslb = brois(1).rois(2).slice;      % Femur bone
   rslf = intersect(rslc,rslb);        % Ensure unique slices in sorted order
%
% Get Tibia Slices
%
   rslc = brois(2).rois(1).slice;      % Tibia cartilage
   rslb = brois(2).rois(2).slice;      % Tibia bone
   rsl = intersect(rslc,rslb);         % Ensure unique slices in sorted order
   rsl = intersect(rsl,rslf);          % Ensure femur and tibia slices
   nrsl = size(rsl,1);
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_2ROIs1.ps'];          % ROI lines print file name
   pnam2 = [fs '_2ROIs2.ps'];          % ROI areas print file name
   pnam3 = [fs '_2ROIs3.ps'];          % ROI regions cylinders print file name
   pnam4 = [fs '_2ROIs4.ps'];          % ROI regions masks print file name
%
   f = cell(2,nrsl);    % Femur coordinates (1 - cartilage, 2 - bone)
   t = cell(2,nrsl);    % Tibia coordinates (1 - cartilage, 2 - bone)
%
   f3 = cell(2,nrsl);   % Femur coordinates (1 - cartilage, 2 - bone)
   t3 = cell(2,nrsl);   % Tibia coordinates (1 - cartilage, 2 - bone)
%
   ibone = false(nrsl,2);              % Bone exists? (1 - femur, 2 - tibia)
   icmprt = ones(nrsl,1);              % Tibal compartments (1 - lateral, 2 - medial)
   npxt = prod(iszs);                  % Total number of pixels in a slice
   maskf = false(npxt,2,nrsl);         % Mask for femoral cartilage
   maskt = false(npxt,2,nrsl);         % Mask for tibial cartilage
%
   cid = cell(nrsl,1);                 % Index to tibia contact points
   cdatt = cell(nrsl,1);               % Tibia contact points
   cdatf = cell(nrsl,1);               % Femur contact points
%
   xyz2 = zeros(2,3);   % Unloaded posterior femoral center points
   vec2 = zeros(2,3);   % Unloaded posterior femoral normal vectors
%
% Loop through Tibia Slices
%
   for k = 1:nrsl       % Loop through slices
%
% Get T2 Slice Image
%
      slk = rsl(k);        % Slice number
      img = squeeze(v(:,:,slk,id5m));
%
      figure;
      orient landscape;
      imagesc(img);
      colormap gray;
      axis image;
      axis off;
      title({dirstr; [fs ' Slice ' int2str(slk)]},'FontSize',16, ...
            'FontWeight','bold');
      hold on;
%
      lh = gobjects(4,1);              % Line graphic handles
      idl = false(4,1);
%
% Get ROI Data for this Slice and Plot ROIs
%
      for lb = 1:2      % Loop through bones - femur = 1 and tibia = 2
         for lc = 1:2   % Loop through surfaces - cartilage = 1 and bone = 2
            idxl = 2*lb+lc-2;
            idxr = brois(lb).rois(lc).slice==slk;
            if any(idxr)
              if lb==1&&itroch
                nc = 3;
              else
                nc = 2;
              end
              for n = 1:nc   % Loop through compartments
                 idxs = brois(lb).rois(lc).roi(n).imageno==slk;
                 if any(idxs)
                   dat = cell2mat( ...
                                 brois(lb).rois(lc).roi(n).data(idxs)');
                   dat3 = cell2mat( ...
                                brois(lb).rois(lc).roi(n).data3(idxs)');
                   [dat,idx] = fix_gap(dat);
                   dat3 = dat3(idx,:);
                   if lb==1
                     f{lc,k} = dat;    % Femur
                     f3{lc,k} = dat3;  % Femur
                   else
                     t{lc,k} = dat;    % Tibia
                     t3{lc,k} = dat3;  % Tibia
                     icmprt(k) = n;
                   end
                   lh(idxl) = plot(dat(:,1),dat(:,2),lt(idxl,:));
                   idl(idxl) = true;
                 end
              end       % End of n loop - lateral/medial loop
            end         % End of if - any rois for this slice
         end            % End of lc loop - cartilage/bone loop
         if lb==1
           ibone(k,lb) = all(~cellfun('isempty',f(:,k)));
         else
           ibone(k,lb) = all(~cellfun('isempty',t(:,k)));
         end
      end               % End of lb loop - femur/tibia loop
%
% Get Transform Parameters
% Scale, s, Rotation Matrix, r, and Translation Vector, tv.
%
      if k==1
        [~,~,r,tv] = trnsf2pixel(f3{2,k},f{2,k},f3{2,k});
        s = 1./pspcs(1);               % Scale
        r3 = -r(:,[1 3 2]);
      end
%
% Get Contact Points
%
      [cid{k},cdatf{k}] = pwl_contct(f{1,k},t{1,k},thresh);
      if ~isempty(cid{k})
        cdatt{k} = t3{1,k}(cid{k},:);
      end
%
% Add Legends and Print Slice Plots
%
      legend(lh(idl),legds(idl,:),'Interpreter','none');
      if k==1
        print('-dpsc2','-r600','-fillpage',pnam1);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam1);
      end
%
% Create Logical Masks for the Cartilage on this Slice
%
      if ibone(k,1)
        [maskf(:,1,k),maskf(:,2,k)] = cr_mask2f(f(:,k),iszs,dist,pspcs);
      end
%
      if ibone(k,2)
        [maskt(:,1,k),maskt(:,2,k)] = cr_mask2(t(:,k),iszs,dist,pspcs);
      end
%
% Plot ROIs
%
      mask1 = maskf(:,1,k)|maskt(:,1,k);    % Superficial cartilage mask
      mask2 = maskf(:,2,k)|maskt(:,2,k);    % Deep cartilage mask
      mask = mask1|mask2;              % All cartilage mask
      cmx = max(img(:));
      img1 = img-cmx-1;
      idmsk = ~mask;
      img2 = img;
      img2(idmsk) = img(idmsk)-max(img(:))-1;
      dcmx = 16*cmx/128;
      img1(mask1) = dcmx;  % Blue - Superficial
      img1(mask2) = cmx-dcmx;   % Red - Deep
%
      figure;
      orient landscape;
      imagesc(img1);
      colormap(cmap);
      caxis([-cmx cmx]);
      axis image;
      axis off;
      title({dirstr; [fs ' Slice ' int2str(slk)]},'FontSize',16, ...
            'FontWeight','bold');
%
      if k==1
        print('-dpsc2','-r600','-fillpage',pnam2);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam2);
      end
%
   end                  % End of k loop - tibia slices
%
% Close Slice Plots
%
   close all;
%
% Get Centers of ROIs from the Tibial Contact Points
%
   ilat = icmprt==1;    % Lateral compartment
   xyzls = cell2mat(cdatt(ilat));      % Lateral tibia XYZ coordinates in mm
   xyzms = cell2mat(cdatt(~ilat));     % Medial tibia XYZ coordinates in mm
%
%    idv = ~cellfun('isempty',cdatt);
%    nt = cellfun('size',cdatt,1);       % Number of contact points on each slice
%
   xyz1(2,:) = mean(xyzms);            % Center of medial contact ROI
   xyz1(1,:) = mean(xyzls);            % Center of lateral contact ROI
%
% Get Centers and Axial Vectors of Posterior Femoral ROIs
% Row 1 -> Lateral, and Row 2 -> Medial
%
   idc2 = fcmprt(fids);
   xyz2(idc2,:) = fpts(fids,:);        % Centers of cylinders
   vec2(idc2,:) = fvec(fids,:);        % Cylinder axes
%
% Transform Cartilage and Femoral Segmentations to ROI Cylinders Axes
%
   ft1 = cell(2,1);     % Loaded ROIs femur-tibia bone transformed (1 - lateral, 2 - medial)
   r2 = zeros(3,3,2);   % Rotation matrix to cylinders axes
   ft2 = cell(2,1);     % Unloaded ROIs femur transformed (1 - lateral, 2 - medial)
%
   for l = 1:2          % Loop through compartments
      ilat = icmprt==l;
%       ft1{l} = coord_tf(xyz1(l,:),r1,[f3(2,ilat); t3(2,ilat)]);
      ft1{l} = coord_tf(xyz1(l,:),r3,[f3(2,ilat); t3(2,ilat)]);
      r2(:,:,l) = vec2rot3d(vec2(l,:),xvec);
      ft2{l} = coord_tf(xyz2(l,:),r2(:,:,l),f3(:,ilat));
   end
%
% Index to Compartment Slices
%
   islc = zeros(nrsl,1);
%
   idcl = find(icmprt==1);             % Lateral compartment
   ncl = size(idcl,1);
   islc(idcl) = (1:ncl)';
%
   idcl = find(icmprt==2);             % Medial compartment
   ncl = size(idcl,1);
   islc(idcl) = (1:ncl)';
%
%  Setup Figures
%
   hf3 = gobjects(5,1);  % Figure handles for rotated cylindrical ROIs
   hf3(5) = figure;         % Figure handle for cylindrical ROIs
   orient landscape;
%
% Plot Cylinders in MRI Coordinates
%
%    xyzct = tf_coord(r1,xyz1(1,:),{xyzc});        % Loaded lateral ROI
   xyzct = tf_coord(r3',xyz1(1,:),{xyzc});       % Loaded lateral ROI
   xyzct = xyzct{1};
   xct = reshape(xyzct(:,1),nrow,ncol);
   yct = reshape(xyzct(:,2),nrow,ncol);
   zct = reshape(xyzct(:,3),nrow,ncol);
%
   plot3(xct,yct,zct,'y-','LineWidth',1);
   hold on;
   plot3(xct',yct',zct','y-','LineWidth',1);
   axis equal;
%
%    xyzct = tf_coord(r1,xyz1(2,:),{xyzc});        % Loaded medial ROI
   xyzct = tf_coord(r3',xyz1(2,:),{xyzc});       % Loaded medial ROI
   xyzct = xyzct{1};
   xct = reshape(xyzct(:,1),nrow,ncol);
   yct = reshape(xyzct(:,2),nrow,ncol);
   zct = reshape(xyzct(:,3),nrow,ncol);
%
   plot3(xct,yct,zct,'y-','LineWidth',1);
   plot3(xct',yct',zct','y-','LineWidth',1);
%
   xyzct = tf_coord(r2(:,:,1)',xyz2(1,:),{xyzc});     % Unloaded lateral ROI
   xyzct = xyzct{1};
   xct = reshape(xyzct(:,1),nrow,ncol);
   yct = reshape(xyzct(:,2),nrow,ncol);
   zct = reshape(xyzct(:,3),nrow,ncol);
%
   plot3(xct,yct,zct,'y-','LineWidth',1);
   plot3(xct',yct',zct','y-','LineWidth',1);
%
   xyzct = tf_coord(r2(:,:,2)',xyz2(2,:),{xyzc});     % Unloaded medial ROI
   xyzct = xyzct{1};
   xct = reshape(xyzct(:,1),nrow,ncol);
   yct = reshape(xyzct(:,2),nrow,ncol);
   zct = reshape(xyzct(:,3),nrow,ncol);
%
   plot3(xct,yct,zct,'y-','LineWidth',1);
   plot3(xct',yct',zct','y-','LineWidth',1);
%
% Plot Cylinders in Rotated Coordinates Aligned with the Axis of the
% Cylinders
%
   hf3(1) = figure;     % Lateral loaded ROI
   orient landscape;
%
   plot3(xc,yc,zc,'y-');
   hold on;
   plot3(xc',yc',zc','y-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Lateral Loaded ROI','FontSize',14,'FontWeight','bold');
%
   hf3(2) = figure;     % Medial loaded ROI
   orient landscape;
%
   plot3(xc,yc,zc,'y-');
   hold on;
   plot3(xc',yc',zc','y-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Medial Loaded ROI','FontSize',14,'FontWeight','bold');
%
   hf3(3) = figure;     % Lateral unloaded ROI
   orient landscape;
%
   plot3(xc,yc,zc,'y-');
   hold on;
   plot3(xc',yc',zc','y-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Lateral Unloaded ROI','FontSize',14,'FontWeight','bold');
%
   hf3(4) = figure;     % Medial unloaded ROI
   orient landscape;
%
   plot3(xc,yc,zc,'y-');
   hold on;
   plot3(xc',yc',zc','y-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Medial Unloaded ROI','FontSize',14,'FontWeight','bold');
%
% Loop through the Slices to Define Analysis Regions
%
   irsl1 = false(nrsl,2);    % Loaded contact ROIs on slice (column 1 - lateral, column 2 - medial)
   irsl2 = false(nrsl,2);    % Unloaded posterior ROIs on slice
%
   maskr1 = false(npxt,2,nrsl);        % Mask for loaded ROI
   maskr2 = false(npxt,2,nrsl);        % Mask for unloaded ROI
%
   for k = 1:nrsl       % Loop through slices
%
% Plot Unrotated Segmentations
%
      figure(hf3(5));
%
      plot3(f3{1,k}(:,1),f3{1,k}(:,2),f3{1,k}(:,3),lt(1,:), ...
            'LineWidth',1);
      plot3(f3{2,k}(:,1),f3{2,k}(:,2),f3{2,k}(:,3),lt(2,:), ...
            'LineWidth',1);
      plot3(t3{1,k}(:,1),t3{1,k}(:,2),t3{1,k}(:,3),lt(3,:), ...
            'LineWidth',1);
      plot3(t3{2,k}(:,1),t3{2,k}(:,2),t3{2,k}(:,3),lt(4,:), ...
            'LineWidth',1);
%
% Find Lateral and Medial Loaded and Unloaded ROIs on Slice
%
      l = icmprt(k);    % Compartments
      ks = islc(k);
%
% Get Coordinates for Both Lines
%
      fb1 = ft1{l}{1,ks};    % Slice segmentations in loaded ROI coordinates
      tb1 = ft1{l}{2,ks};
%
      fb2 = ft2{l}{1,ks};    % Slice segmentations in unloaded ROI coordinates
      tb2 = ft2{l}{2,ks};
%
% Get Polygon ROIs in Rotated Coordinates Based on Intersections with
% Cylindrical ROIs
%
      [xyzp1,xyzif1,xyzit1] = cyl_pwl(fb1,tb1,rrad,15,24);
      [xyzp2,xyzif2,xyzit2] = cyl_pwl(fb2,tb2,rrad,15,12);
%
% Intersections with the Loaded Cylindrical ROIs?
%
      if ~isempty(xyzp1)
%
        irsl1(k,l) = true;   % Intersection true for this slice and compartment
%
% Plot Lines and Polygon
%
        figure(hf3(l));
%
        plot3(fb1(:,1),fb1(:,2),fb1(:,3),'b.-');
        plot3(tb1(:,1),tb1(:,2),tb1(:,3),'g.-');
%
        plot3(fb1(1,1),fb1(1,2),fb1(1,3),'b^');            % Start
        plot3(fb1(end,1),fb1(end,2),fb1(end,3),'bs');      % End
%
        plot3(tb1(1,1),tb1(1,2),tb1(1,3),'g^');            % Start
        plot3(tb1(end,1),tb1(end,2),tb1(end,3),'gs');      % End
%
        if ~isempty(xyzif1)
          plot3(xyzif1(:,1),xyzif1(:,2),xyzif1(:,3),'rs');
        end
        if ~isempty(xyzit1)
          plot3(xyzit1(:,1),xyzit1(:,2),xyzit1(:,3),'rs');
        end
%
        hp = patch(xyzp1(:,1),xyzp1(:,2),xyzp1(:,3),'k');
        set(hp,'FaceAlpha',0.4);
        set(hp,'EdgeColor','none');
%
% Transform Back to MRI Coordinates
%
%         xyzpt1 = tf_coord(r1,xyz1(l,:),{xyzp1});
        xyzpt1 = tf_coord(r3',xyz1(l,:),{xyzp1});
        xyzpt1 = xyzpt1{1};
%
        figure(hf3(5));
        hp = patch(xyzpt1(:,1),xyzpt1(:,2),xyzpt1(:,3),'k');
        set(hp,'FaceAlpha',0.4);
        set(hp,'EdgeColor','none');
if l==2
keyboard
end
%
% Transform to Pixel Coordinates
%
        xyzpp1 = s*xyzpt1*r+repmat(tv,size(xyzpt1,1),1);
%
% Create Mask from Polygon
%
        msk = poly2mask(xyzpp1(:,1),xyzpp1(:,2),iszs(1),iszs(2));
        maskr1(:,l,k) = msk(:);
%
      end
%
% Intersections with the Unloaded Cylindrical ROIs?
%
      if ~isempty(xyzp2)
%
        irsl2(k,l) = true;   % Intersection true for this slice and compartment
%
% Plot Lines and Polygon
%
        figure(hf3(l+2));
%
        plot3(fb2(:,1),fb2(:,2),fb2(:,3),'b.-');
        plot3(tb2(:,1),tb2(:,2),tb2(:,3),'g.-');
%
        plot3(fb2(1,1),fb2(1,2),fb2(1,3),'b^');            % Start
        plot3(fb2(end,1),fb2(end,2),fb2(end,3),'bs');      % End
%
        plot3(tb2(1,1),tb2(1,2),tb2(1,3),'g^');            % Start
        plot3(tb2(end,1),tb2(end,2),tb2(end,3),'gs');      % End
%
        if ~isempty(xyzif2)
          plot3(xyzif2(:,1),xyzif2(:,2),xyzif2(:,3),'rs');
        end
        if ~isempty(xyzit2)
          plot3(xyzit2(:,1),xyzit2(:,2),xyzit2(:,3),'rs');
        end
%
        hp = patch(xyzp2(:,1),xyzp2(:,2),xyzp2(:,3),'k');
        set(hp,'FaceAlpha',0.4);
        set(hp,'EdgeColor','none');
%
% Transform Back to MRI Coordinates
%
        xyzpt2 = tf_coord(r2(:,:,l)',xyz2(l,:),{xyzp2});
        xyzpt2 = xyzpt2{1};
%
        figure(hf3(5));
        hp = patch(xyzpt2(:,1),xyzpt2(:,2),xyzpt2(:,3),'k');
        set(hp,'FaceAlpha',0.4);
        set(hp,'EdgeColor','none');
%
% Transform to Pixel Coordinates
%
        xyzpp2 = s*xyzpt2*r+repmat(tv,size(xyzpt2,1),1);
%
% Create Mask from Polygon
%
        msk = poly2mask(xyzpp2(:,1),xyzpp2(:,2),iszs(1),iszs(2));
        maskr2(:,l,k) = msk(:);
%
      end
%
   end                  % End of slice loop - k
%
% Print Plots
%
   figure(hf3(5));
   title({dirstr; fs; 'T2* ROIs'},'FontSize',12,'FontWeight','bold');
%
   for k = 1:5
%
      figure(hf3(k));
%
      if k==1
        print('-dpsc2','-r600','-fillpage',pnam3);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam3);
      end
%
      if k==5
        view(2);
        print('-dpsc2','-r600','-fillpage','-append',pnam3);
      end
%
   end
%
   close all;
%
% Combine Cartilage and ROI Masks to Get Lateral Loaded Mask
%
   rsll1 = rsl(irsl1(:,1));
   nrsll1 = size(rsll1,1);
%
   maskfrl1 = false(npxt,2,nrsll1);    % Mask for loaded region femoral cartilage
   masktrl1 = false(npxt,2,nrsll1);    % Mask for loaded region tibial cartilage
%
   for k = 1:nrsll1
%
      ida = ismember(rsl,rsll1(k));
%
      for l = 1:2       % 1 - Superficial, 2 - Deep
%
         maskfrl1(:,l,k) = maskf(:,l,ida)&maskr1(:,1,ida);
         masktrl1(:,l,k) = maskt(:,l,ida)&maskr1(:,1,ida);
%
      end
%
   end
%
% Combine Cartilage and ROI Masks to Get Medial Loaded Mask
%
   rslm1 = rsl(irsl1(:,2));
   nrslm1 = size(rslm1,1);
%
   maskfrm1 = false(npxt,2,nrslm1);    % Mask for loaded region femoral cartilage
   masktrm1 = false(npxt,2,nrslm1);    % Mask for loaded region tibial cartilage
%
   for k = 1:nrslm1
%
      ida = ismember(rsl,rslm1(k));
%
      for l = 1:2       % 1 - Superficial, 2 - Deep
%
         maskfrm1(:,l,k) = maskf(:,l,ida)&maskr1(:,2,ida);
         masktrm1(:,l,k) = maskt(:,l,ida)&maskr1(:,2,ida);
%
      end
%
   end
%
% Combine Cartilage and ROI Masks to Get Lateral Unloaded Femoral Mask
%
   rsll2 = rsl(irsl2(:,1));
   nrsll2 = size(rsll2,1);
%
   maskfrl2 = false(npxt,2,nrsll2);    % Mask for unloaded region femoral cartilage
%
   for k = 1:nrsll2
%
      ida = ismember(rsl,rsll2(k));
%
      for l = 1:2       % 1 - Superficial, 2 - Deep
%
         maskfrl2(:,l,k) = maskf(:,l,ida)&maskr2(:,1,ida);
%
      end
%
   end
%
% Combine Cartilage and ROI Masks to Get Medial Unloaded Femoral Mask
%
   rslm2 = rsl(irsl2(:,2));
   nrslm2 = size(rslm2,1);
%
   maskfrm2 = false(npxt,2,nrslm2);    % Mask for unloaded region femoral cartilage
%
   for k = 1:nrslm2
%
      ida = ismember(rsl,rslm2(k));
%
      for l = 1:2       % 1 - Superficial, 2 - Deep
%
         maskfrm2(:,l,k) = maskf(:,l,ida)&maskr2(:,2,ida);
%
      end
%
   end
%
% Plot the ROI Masks
%
   rslp = unique([rsll1; rsll2; rslm1; rslm2]);
   nrslp = size(rslp,1);
%
   for k = 1:nrslp
%
      slk = rslp(k);
      idc = rsl==slk;
      idc = icmprt(idc);
%
      img = squeeze(v(:,:,slk,id5m));
      cmx = max(img(:));
      img1 = img-cmx-1;
%
      figure;
      orient landscape;
%
% Plot Loaded Contact and Unloaded Posterior ROIs
%
      mask1 = false(npxt,1);
      mask2 = false(npxt,1);
%
      idp = rsll1==slk;
      if any(idp)
        mask1 = maskfrl1(:,1,idp)|masktrl1(:,1,idp)|mask1; % Superficial cartilage mask
        mask2 = maskfrl1(:,2,idp)|masktrl1(:,2,idp)|mask2; % Deep cartilage mask
      end
%
      idp = rslm1==slk;
      if any(idp)
        mask1 = maskfrm1(:,1,idp)|masktrm1(:,1,idp)|mask1; % Superficial cartilage mask
        mask2 = maskfrm1(:,2,idp)|masktrm1(:,2,idp)|mask2; % Deep cartilage mask
      end
%
      idp = rsll2==slk;
      if any(idp)
        mask1 = maskfrl2(:,1,idp)|mask1;    % Superficial cartilage mask
        mask2 = maskfrl2(:,2,idp)|mask2;    % Deep cartilage mask
      end
%
      idp = rslm2==slk;
      if any(idp)
        mask1 = maskfrm2(:,1,idp)|mask1;    % Superficial cartilage mask
        mask2 = maskfrm2(:,2,idp)|mask2;    % Deep cartilage mask
      end
%
      dcmx = 16*cmx/128;
      img1(mask1) = dcmx;         % Blue - Superficial
      img1(mask2) = cmx-dcmx;     % Red - Deep
%
      imagesc(img1);
      colormap(cmap);
      caxis([-cmx cmx]);
      axis image;
      axis off;
      title({[fs ' Slice ' int2str(slk)]; [cmprt{idc}, ...
            ' Compartment']; 'T2* ROIs'}, ...
            'FontSize',12,'FontWeight','bold');
%
      if k==1
        print('-dpsc2','-r600','-fillpage',pnam4);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam4);
      end
%
   end
%
% Save Masks, ROIs and Slice Information into MAT File
%
   savnam = [mnam(1:end-4) '_2rois.mat'];
   save(savnam,'f','f3','ibone','icmprt','maskf','maskfrl1', ...
               'maskfrm1','maskfrl2','maskfrm2','maskr1','maskr2', ...
               'maskt','masktrl1','masktrm1','brois','rsl', ...
               'rsll1','rslm1','rsll2','rslm2','t','t3');
%
   close all;
%
end                     % End of m loop - MAT file loop
%
return