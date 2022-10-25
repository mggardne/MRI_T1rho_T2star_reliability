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
%             3.  The M-files cr_mask2.m, cr_mask2f.m, cyl_lin3.m,
%             cyl_pwl.m, fix_gap.m, in_tri2d.m,
%             lsect2.m, lsect2a.m, lsect3.m, lsect4.m, lsect5.m,
%             midline.m, mk2_tri_2d.m, mk2_tri_2df.m, pts2lin.m,
%             mk2_tri_2d.m, pts2lin.m, pwl_contct.m, rd_rois3.m and
%             rd_roi6.m must be in the current directory or path.
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
% Axis Vector and Contact ROIs Rotation Matrix
%
xvec = [1 0 0];         % X-vector for defining 3D cylinder coordinate systems
r1 = eye(3);            % No rotation
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
% T1rho
%
rdir = 'RHO';           % Directory for T1rho segmentations
id5m = 1;               % Use first spin lock time for plots
thresh = 1.953125;      % Four (4) pixels = 4*0.48828125 = 1.953125
% threshs = num2str(thresh,'%.3f');      % Threshold as a string
fr = find(fres==0);
fids = intersect(fidx,fr);             % Subject, visit and T1rho analysis
%
% Get T1rho MAT Files in Directory
%
d = dir('T1rho_S*.mat');
mnams = {d.name}';
idr = contains(mnams,'roi','IgnoreCase',true);
mnams = mnams(~idr);
nmat = size(mnams,1);
%
% Loop through T1rho MAT Files
%
% for m = 1:nmat
for m = 1:1
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
   fids = intersect(fids,fl);          % Include index to leg
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
   rslmn = min(rslb);   % Tibia bone minimum slice
   rslmx = max(rslb);   % Tibia bone maximum slice
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_2ROIs1.ps'];          % ROI lines print file name
   pnam2 = [fs '_2ROIs2.ps'];          % ROI areas print file name
   pnam3 = [fs '_2ROIs3.ps'];          % ROI regions cylinders print file name
   pnam4 = [fs '_2ROIs3.ps'];          % ROI regions masks print file name
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
   vec1 = repmat([0 0 1],2,1);         % Cylinder axis is inferior-superior in MRI coordinate system
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
      ft1{l} = coord_tf(xyz1(l,:),r1,[f3(2,ilat); t3(2,ilat)]);
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
%  Setup Figure
%
   hf3 = figure;
   orient landscape;
%
% Plot Cylinders
%
   [xc,yc,zc] = cylinder(rrad*ones(4,1),72);
   zc = 45*zc-15;       % Scale and translate Z coordinates
%
   hs = gobjects(4,1);  % Subplot handles
%
   hs(1) = subplot(2,2,1);   % Lateral loaded ROI
   plot3(xc,yc,zc,'k-');
   hold on;
   plot3(xc',yc',zc','k-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Lateral Loaded ROI','FontSize',14,'FontWeight','bold');
%
   hs(2) = subplot(2,2,2);   % Medial loaded ROI
   plot3(xc,yc,zc,'k-');
   hold on;
   plot3(xc',yc',zc','k-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Medial Loaded ROI','FontSize',14,'FontWeight','bold');
%
   zc = zc-15;          % Translate Z coordinates
   hs(3) = subplot(2,2,3);   % Lateral unloaded ROI
   plot3(xc,yc,zc,'k-');
   hold on;
   plot3(xc',yc',zc','k-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Lateral Unloaded ROI','FontSize',14,'FontWeight','bold');
%
   hs(4) = subplot(2,2,4);   % Medial unloaded ROI
   plot3(xc,yc,zc,'k-');
   hold on;
   plot3(xc',yc',zc','k-');
   axis equal;
   zlabel('Z (mm)','FontSize',11,'FontWeight','bold');
   title('Medial Unloaded ROI','FontSize',14,'FontWeight','bold');
%
% Loop through the Slices to Define Analysis Regions
%
   irsl1 = false(nrsl,2);    % Contact ROI on slice (column 1 - lateral, column 2 - medial)
   irsl2 = false(nrsl,2);    % Posterior ROI on slice
%
   maskr1 = false(npxt,2,nrsl);        % Mask for loaded ROI
   maskr2 = false(npxt,2,nrsl);        % Mask for unloaded ROI
%
   for k = 1:nrsl       % Loop through slices
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
      [xyzp1,xyzif1,xyzit1] = cyl_pwl(fb1,tb1,rrad);
      [xyzp2,xyzif2,xyzit2] = cyl_pwl(fb2,tb2,rrad);
%
% Intersections with the Loaded Cylindrical ROIs?
%
      if ~isempty(xyzp1)
%
        irsl1(k,l) = true;   % Intersection true for this slice and compartment
%
% Plot Lines and Polygon
%
        subplot(hs(l));
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
        plot3(xyzif1(:,1),xyzif1(:,2),xyzif1(:,3),'rs');
        plot3(xyzit1(:,1),xyzit1(:,2),xyzit1(:,3),'rs');
%
        hp = patch(xyzp1(:,1),xyzp1(:,2),xyzp1(:,3),'k');
        set(hp,'FaceAlpha',0.5);
        set(hp,'EdgeColor','none');
%
% Transform Back to MRI Coordinates
%
        xyzpt1 = tf_coord(r1,xyz1(l,:),{xyzp1});
        xyzpt1 = xyzpt1{1};
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
        subplot(hs(l+2));
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
        plot3(xyzif2(:,1),xyzif2(:,2),xyzif2(:,3),'rs');
        plot3(xyzit2(:,1),xyzit2(:,2),xyzit2(:,3),'rs');
%
        hp = patch(xyzp2(:,1),xyzp2(:,2),xyzp2(:,3),'k');
        set(hp,'FaceAlpha',0.5);
        set(hp,'EdgeColor','none');
%
% Transform Back to MRI Coordinates
%
        xyzpt2 = tf_coord(r2(:,:,l),xyz2(l,:),{xyzp2});
        xyzpt2 = xyzpt2{1};
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
   print('-dpsc2','-r600','-fillpage',pnam3);
%
   for k = 1:4
      subplot(hs(k));
      view(-90,90);
   end
%
   print('-dpsc2','-r600','-fillpage','-append',pnam3);
%
keyboard
%
% Get Lateral ROI Masks
%
   slliml = fix(slliml./pspcs(1));
   maskfrl = false(npxt,2,nrsll);      % Mask for lateral femoral cartilage
   masktrl = false(npxt,2,nrsll);      % Mask for lateral tibial cartilage
%
   for k = 1:nrsll
%
      slk = rsll(k);
      maskrl = false(iszs);
      maskrl(:,slliml(k,1):slliml(:,2)) = true;
      maskrl = repmat(maskrl(:),1,2);
      idl = slk==rsl;
      maskfrl(:,:,k) = maskf(:,:,idl)&maskrl;
      masktrl(:,:,k) = maskt(:,:,idl)&maskrl;
%
   end
%
% Get Medial ROI Masks
%
   sllimm = fix(sllimm./pspcs(1));
   maskfrm = false(npxt,2,nrslm);      % Mask for medial femoral cartilage
   masktrm = false(npxt,2,nrslm);      % Mask for medial tibial cartilage
%
   for k = 1:nrslm
%
      slk = rslm(k);
      maskrm = false(iszs);
      maskrm(:,sllimm(k,1):sllimm(:,2)) = true;
      maskrm = repmat(maskrm(:),1,2);
      idm = slk==rsl;
      maskfrm(:,:,k) = maskf(:,:,idm)&maskrm;
      masktrm(:,:,k) = maskt(:,:,idm)&maskrm;
%
   end
%
% Plot the ROI Masks
%
   close all;
%
   rslp = unique([rsll; rslm]);
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
      orient tall;
      subplot(2,1,1);
%
% Plot Contact ROIs on the Top Half of the Page
%
      if any(rsll==slk)||any(rslm==slk)
        if any(rsll==slk)
          idp = rsll==slk;
          mask1 = maskfrl(:,1,idp)|masktrl(:,1,idp);  % Superficial cartilage mask
          mask2 = maskfrl(:,2,idp)|masktrl(:,2,idp);  % Deep cartilage mask
          dcmx = 16*cmx/128;
          img1(mask1) = dcmx;  % Blue - Superficial
          img1(mask2) = cmx-dcmx;   % Red - Deep
        else
          idp = rslm==slk;
          mask1 = maskfrm(:,1,idp)|masktrm(:,1,idp);  % Superficial cartilage mask
          mask2 = maskfrm(:,2,idp)|masktrm(:,2,idp);  % Deep cartilage mask
          dcmx = 16*cmx/128;
          img1(mask1) = dcmx;  % Blue - Superficial
          img1(mask2) = cmx-dcmx;   % Red - Deep
        end
      end
%
      imagesc(img1);
      colormap(cmap);
      caxis([-cmx cmx]);
      axis image;
      axis off;
      title({[fs ' Slice ' int2str(slk)]; [cmprt{idc}, ...
            ' Compartment']; 'T1\rho Contact Points ROIs'}, ...
            'FontSize',12,'FontWeight','bold');
%
      print('-dpsc2','-r600','-fillpage','-append',pnam4);
%
   end
%
% Save Masks, ROIS and Slice Information into MAT File
% Note maskp and p are empty/false.
%
   savnam = [mnam(1:end-4) '_2rois.mat'];
   save(savnam,'f','ibone','icmprt','maskf', ...
               'maskfrl','maskfrm','maskt', ...
               'masktrl','masktrm','brois', ...
               'rsl','rsll','rslm','t');
%
   close all;
%
end                     % End of m loop - MAT file loop
%
return
%
% T2star
%
rdir = 'T2S';           % Directory for T2* segmentations
id5 = repmat(3,4,1);    % Default echo time for segmentations
thresh = 1.75;          % Four (4) pixels = 4*0.4375 = 1.75
% threshs = num2str(thresh,'%.3f');      % Threshold as a string
fr = find(fres==1);
fids = intersect(fidx,fr);             % Subject, visit and T2* analysis
%
% Get T2* MAT Files in Directory
%
d = dir('T2star_S*.mat');
mnams = {d.name}';
idr = contains(mnams,'roi','IgnoreCase',true);
mnams = mnams(~idr);
nmat = size(mnams,1);
%
% Loop through T2* MAT Files
%
% for m = 1:nmat
for m = 1:1
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
   fids = intersect(fids,fl);          % Include index to leg
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
   brois = rd_rois(rdir,leg,ld,itroch,5);
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
   rsl3 = 2.0*rsl;      % Slice thickness = 2 mm
%
   rslmn = min(rslb);   % Tibia bone minimum slice
   rslmx = max(rslb);   % Tibia bone maximum slice
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_2ROIs1.ps'];          % ROI lines print file name
   pnam2 = [fs '_2ROIs2.ps'];          % ROI areas print file name
   pnam3 = [fs '_2ROIs3.ps'];          % ROI regions print file name
%
   f = cell(2,nrsl);    % Femur coordinates (1 - cartilage, 2 - bone)
   t = cell(2,nrsl);    % Tibia coordinates (1 - cartilage, 2 - bone)
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
   ys = cell(nrsl,1);   % Y coordinates for contact points
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
              for n = 1:nc    % Loop through compartments
                 idxs = brois(lb).rois(lc).roi(n).imageno==slk;
                 if any(idxs)
                   dat = cell2mat(brois(lb).rois(lc).roi(n).data(idxs)');
                   dat = fix_gap(dat);
                   if lb==1
                     f{lc,k} = dat;  % Femur
                   else
                     t{lc,k} = dat;  % Tibia
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
% Get Contact Points
%
      [cid{k},cdatf{k}] = pwl_contct(f{1,k},t{1,k},thresh);
      if ~isempty(cid{k})
        cdatt{k} = t{1,k}(cid{k},:);
        ys{k} = repmat(rsl3(k),size(cdatt{k},1),1);
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
% Get Centers of ROIs from the Contact Points
%
   ilat = icmprt==1;    % Lateral compartment
   xzl = cell2mat(cdatt(ilat))*pspcs(1);    % Lateral XZ coordinates - convert to mm
   xzm = cell2mat(cdatt(~ilat))*pspcs(1);   % Medial XZ coordinates - convert to mm
%
%    idv = ~cellfun('isempty',cdatt);
%    nt = cellfun('size',cdatt,1);       % Number of contact points on each slice
%
   yl = cell2mat(ys(ilat));            % Lateral Y coordinates
   ym = cell2mat(ys(~ilat));           % Medial Y coordinates
%
   xyzls = [xzl(:,1) yl xzl(:,2)];     % 3D lateral contact points
   xyzl = mean(xyzls);                 % Center of lateral contact ROI
   xyzms = [xzm(:,1) ym xzm(:,2)];     % 3D medial contact points
   xyzm = mean(xyzms);                 % Center of medial contact ROI
%
% Get Slices Within RRAD of the Contact Points ROIs
%
   idxl = (rsl3>xyzl(2)-rrad)&(rsl3<xyzl(2)+rrad);    % Lateral
   rsll = rsl(idxl);    % Lateral compartment slices
   idxm = (rsl3>xyzm(2)-rrad)&(rsl3<xyzm(2)+rrad);    % Medial
   rslm = rsl(idxm);    % Medial compartment slices
%
%  Setup Figure
%
   hf3 = figure;
   orient landscape;
   hold on;
%
% Loop through the Slices to Define Analysis Regions
%
   idl = 0;
   idm = 0;
   ifirst1 = true;
   ifirst2 = true;
%
   nrsll = size(rsll,1);     % Number of lateral slices in contact ROI
   nrslm = size(rslm,1);     % Number of medial slices in contact ROI
%
   slliml = zeros(nrsll,2);  % Column 1 - minimum, column 2 - maximum
   sllimm = zeros(nrslm,2);  % Column 1 - minimum, column 2 - maximum
%
   for k = 1:nrsl       % Loop through slices
      slk = rsl(k);
      if ibone(k,2)
        xyzc = t{1,k}*pspcs(1);        % Assumes square pixels
        npts = size(xyzc,1);
        xyzc = [xyzc(:,1) repmat(2.0*slk,npts,1) xyzc(:,2)];
        xyzb = t{2,k}*pspcs(1);        % Assumes square pixels
        npts = size(xyzb,1);
        xyzb = [xyzb(:,1) repmat(2.0*slk,npts,1) xyzb(:,2)];
        figure(hf3);
        plot3(xyzb(:,1),xyzb(:,2),xyzb(:,3),lt(4,:),'LineWidth',1);
        plot3(xyzc(:,1),xyzc(:,2),xyzc(:,3),lt(3,:),'LineWidth',1);
%
        if any(rsll==slk)||any(rslm==slk)
%
          xyzcf = f{1,k}*pspcs(1);     % Assumes square pixels
          npts = size(xyzcf,1);
          if npts>0
            xyzcf = [xyzcf(:,1) repmat(1.5*slk,npts,1) xyzcf(:,2)];
            plot3(xyzcf(:,1),xyzcf(:,2),xyzcf(:,3),lt(1,:), ...
                  'LineWidth',1);
          end
%
          xyzbf = f{2,k}*pspcs(1);     % Assumes square pixels
          npts = size(xyzbf,1);
          xyzbf = [xyzbf(:,1) repmat(1.5*slk,npts,1) xyzbf(:,2)];
          plot3(xyzbf(:,1),xyzbf(:,2),xyzbf(:,3),lt(2,:), ...
                'LineWidth',1);
%
          if any(rsll==slk)            % Lateral
%
            plot3(xyzl(:,1),xyzl(:,2),xyzl(:,3),'ks');     % ROI center
            plot3(xyzls(:,1),xyzls(:,2),xyzls(:,3),'co');  % Contact points
            [xp,yp,zp] = cylinder(repmat(rrad,4,1),72);
            xp1 = xp+xyzl(1);
            yp1 = yp+xyzl(2);
            zp1 = -20.0*zp+xyzl(3)+5.0;
            plot3(xp1,yp1,zp1,'b-','Color',[0 0 0.75],'LineWidth',1);
            plot3(xp1',yp1',zp1','b-','Color',[0 0 0.75],'LineWidth',1);
%
            dy = xyzb(1,2)-xyzl(:,2);
            dy2 = dy*dy;
            dx = sqrt(rrad2-dy2);
            idl = idl+1;
            slliml(idl,:) = [xyzl(:,1)-dx xyzl(:,1)+dx];
%
          else          % Medial
%
            plot3(xyzm(:,1),xyzm(:,2),xyzm(:,3),'ks');     % ROI center
            plot3(xyzms(:,1),xyzms(:,2),xyzms(:,3),'co');  % Contact points
            [xp,yp,zp] = cylinder(repmat(rrad,4,1),72);
            xp1 = xp+xyzm(1);
            yp1 = yp+xyzm(2);
            zp1 = -20.0*zp+xyzm(3)+5.0;
            plot3(xp1,yp1,zp1,'g-','Color',[0 0.6 0],'LineWidth',1);
            plot3(xp1',yp1',zp1','g-','Color',[0 0.6 0],'LineWidth',1);
%
            dy = xyzb(1,2)-xyzm(:,2);
            dy2 = dy*dy;
            dx = sqrt(rrad2-dy2);
            idm = idm+1;
            sllimm(idm,:) = [xyzm(:,1)-dx xyzm(:,1)+dx];
%
          end
        end             % End of if slice in lateral or medial region
%
      end               % End of if tibia bone
   end                  % End of slice loop - k
%
   set(gca,'ZDir','reverse');
   view(3);
   axis equal;
   title({dirstr; [fs ' - T2*']; 'Blue - Lateral/Green - Medial'}, ...
          'FontSize',16,'FontWeight','bold');
   print('-dpsc2','-r600','-fillpage',pnam3);
   view(-90,90);
   print('-dpsc2','-r600','-fillpage','-append',pnam3);
%
% Get Lateral ROI Masks
%
   slliml = fix(slliml./pspcs(1));
   maskfrl = false(npxt,2,nrsll);      % Mask for lateral femoral cartilage
   masktrl = false(npxt,2,nrsll);      % Mask for lateral tibial cartilage
%
   for k = 1:nrsll
%
      slk = rsll(k);
      maskrl = false(iszs);
      maskrl(:,slliml(k,1):slliml(:,2)) = true;
      maskrl = repmat(maskrl(:),1,2);
      idl = slk==rsl;
      maskfrl(:,:,k) = maskf(:,:,idl)&maskrl;
      masktrl(:,:,k) = maskt(:,:,idl)&maskrl;
%
   end
%
% Get Medial ROI Masks
%
   sllimm = fix(sllimm./pspcs(1));
   maskfrm = false(npxt,2,nrslm);      % Mask for medial femoral cartilage
   masktrm = false(npxt,2,nrslm);      % Mask for medial tibial cartilage
%
   for k = 1:nrslm
%
      slk = rslm(k);
      maskrm = false(iszs);
      maskrm(:,sllimm(k,1):sllimm(:,2)) = true;
      maskrm = repmat(maskrm(:),1,2);
      idm = slk==rsl;
      maskfrm(:,:,k) = maskf(:,:,idm)&maskrm;
      masktrm(:,:,k) = maskt(:,:,idm)&maskrm;
%
   end
%
% Plot the ROI Masks
%
   close all;
%
   rslp = unique([rsll; rslm]);
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
      orient tall;
      subplot(2,1,1);
%
% Plot Contact ROIs on the Top Half of the Page
%
      if any(rsll==slk)||any(rslm==slk)
        if any(rsll==slk)
          idp = rsll==slk;
          mask1 = maskfrl(:,1,idp)|masktrl(:,1,idp);  % Superficial cartilage mask
          mask2 = maskfrl(:,2,idp)|masktrl(:,2,idp);  % Deep cartilage mask
          dcmx = 16*cmx/128;
          img1(mask1) = dcmx;  % Blue - Superficial
          img1(mask2) = cmx-dcmx;   % Red - Deep
        else
          idp = rslm==slk;
          mask1 = maskfrm(:,1,idp)|masktrm(:,1,idp);  % Superficial cartilage mask
          mask2 = maskfrm(:,2,idp)|masktrm(:,2,idp);  % Deep cartilage mask
          dcmx = 16*cmx/128;
          img1(mask1) = dcmx;  % Blue - Superficial
          img1(mask2) = cmx-dcmx;   % Red - Deep
        end
      end
%
      imagesc(img1);
      colormap(cmap);
      caxis([-cmx cmx]);
      axis image;
      axis off;
      title({[fs ' Slice ' int2str(slk)]; [cmprt{idc}, ...
            ' Compartment']; 'T2* Contact Points ROIs'}, ...
            'FontSize',12,'FontWeight','bold');
%
      print('-dpsc2','-r600','-fillpage','-append',pnam3);
%
   end
%
% Save Masks, ROIs and Slice Information into MAT File
% Note maskp and p are empty/false.
%
   savnam = [mnam(1:end-4) '_2rois.mat'];
   save(savnam,'f','ibone','icmprt','maskf', ...
               'maskfrl','maskfrm','maskt', ...
               'masktrl','masktrm','brois', ...
               'rsl','rsll','rslm','t');
%
   close all;
%
end                     % End of m loop - MAT file loop
%
return