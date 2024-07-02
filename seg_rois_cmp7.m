%#######################################################################
%
%  * SEGmentation to Regions of Interest (ROIs) CoMParison 7 Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     CSV files to create masks for cylindrical regions of interest in
%     the lateral and medial tibial compartments.  The masks are saved
%     in MAT files with the series number and ending in "_rois7.mat."
%     The center of the cylindrical ROIs is determined using two
%     different methods.  One method is based on fixed positions as a
%     percent of the bony tibia surface and the other method is the
%     center of the contact area on the tibia and femur cartilage
%     surface.  Only for T1rho for subject 7 on visit 2.
%
%     NOTES:  1.  The registered MRI MAT files must be in subject
%             directories starting with "MRIR" and either "Visit1" or
%             "Visit2" visit subdirectories.
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  See rd_dicom.m.
%
%             3.  M-file cr_mask2.m, cr_mask2f.m, fix_gap.m, in_tri2d.m,
%             lsect2.m, lsect2a.m, lsect3.m, lsect4.m, lsect5.m,
%             midline.m, mk2_tri_2d.m, mk2_tri_2df.m, pts2lin.m,
%             mk2_tri_2d.m, pts2lin.m, pwl_contct.m, rd_rois7.m and
%             rd_roi6.m must be in the current directory or path.
%
%             4.  Only for T1rho for subject 7 on visit 2.
%
%     26-May-2022 * Mack Gardner-Morse
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
% iskip = true;
iskip = false;
if ~iskip
%
% T1rho
%
id5m = 1;               % Use first spin lock time for plots
rdir = 'RHO';           % Directory for T1rho segmentations
thresh = 1.953125;      % Four (4) pixels = 4*0.48828125 = 1.953125
% threshs = num2str(thresh,'%.3f');      % Threshold as a string
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
for m = 1:nmat
% for m = 4:nmat
% for m = 1
%
   mnam = mnams{m};
   load(mnam);
   fs = ['S' snt];      % Series number prefaced with a 'S'
%
% Parse Series Text for Leg and Load
%
   if strcmpi(st(1),'L')
     leg = 'L';
   else
     leg = 'R';
   end
%
   if contains(st,'Load','IgnoreCase',true)
     ld = 'LD';
   else
     ld = 'UL';
   end
%
% Read ROIs
%
   if  strcmpi(leg,'L')
     brois = rd_rois(rdir,leg,ld,itroch,4);
   else
     brois = rd_rois7(rdir,leg,ld,itroch,4);
   end
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
% Get Slices Within RRAD of the 20% and 80% Width of Tibia Bone
%
   rslctrs = rslmn+(rslmx-rslmn).*[0.2 0.8];     % 1st column = 20%, 2nd column = 80%
   rslc3 = rslctrs*1.5; % Slice thickness = 1.5 mm
   rsl3 = 1.5*rsl;      % Slice thickness = 1.5 mm
%
   idx = (rsl3>rslc3-rrad)&(rsl3<rslc3+rrad);    % 1st column = 20%, 2nd column = 80%
   rsl1 = rsl(idx(:,1));               % First compartment slices at 20%
   rsl2 = rsl(idx(:,2));               % Second compartment slices at 80%
   regpx = zeros(2);    % Minimum and maximum pixel coordinates for both lateral/medial regions
   regpx(1,:) = iszs(1);     % Minimums
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_ROIs1.ps'];           % ROI lines print file name
   pnam2 = [fs '_ROIs2.ps'];           % ROI areas print file name
   pnam3 = [fs '_ROIs3.ps'];           % ROI regions print file name
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
              for n = 1:nc    % Loop through compartments - lateral =1, medial = 2 and trochlea = 3
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
% Check for Region Slices
%
      if any(rsl1==slk)||any(rsl2==slk)
        if any(rsl1==slk)
          icol = 1;
        else
          icol = 2;
        end
        dat2 = t{2,k}(:,1);
        xmn = min(dat2);
        xmx = max(dat2);
        if xmn<regpx(1,icol)
          regpx(1,icol) = xmn;
        end
        if xmx>regpx(2,icol)
          regpx(2,icol) = xmx;
        end
      end
%
   end                  % End of k loop - tibia slices
%
% Close Slice Plots
%
   close all;
%
% Get Centers of Regions
%
   rctrx = regpx(1,:)'+diff(regpx)'/2;
   rctrx3 = rctrx*pspcs(1);
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
% Get Compartment for Each Center
%
   ksl = rsl==rsl1(1);
   irsl1 = icmprt(ksl);
   ksl = rsl==rsl2(end);
   irsl2 = icmprt(ksl);
%
% Loop through the Slices to Define Analysis Regions
%
   id1 = 0;
   id2 = 0;
   idl = 0;
   idm = 0;
   ifirst1 = true;
   ifirst2 = true;
%
   nrsl1 = size(rsl1,1);
   nrsl2 = size(rsl2,1);
   nrsll = size(rsll,1);     % Number of lateral slices in contact ROI
   nrslm = size(rslm,1);     % Number of medial slices in contact ROI
%
   sllim1 = zeros(nrsl1,2);  % Column 1 - minimum, column 2 - maximum
   sllim2 = zeros(nrsl2,2);  % Column 1 - minimum, column 2 - maximum
   slliml = zeros(nrsll,2);  % Column 1 - minimum, column 2 - maximum
   sllimm = zeros(nrslm,2);  % Column 1 - minimum, column 2 - maximum
%
   for k = 1:nrsl       % Loop through slices
      slk = rsl(k);
      if ibone(k,2)
        xyzc = t{1,k}*pspcs(1);        % Assumes square pixels
        npts = size(xyzc,1);
        xyzc = [xyzc(:,1) repmat(1.5*slk,npts,1) xyzc(:,2)];
        xyzb = t{2,k}*pspcs(1);        % Assumes square pixels
        npts = size(xyzb,1);
        xyzb = [xyzb(:,1) repmat(1.5*slk,npts,1) xyzb(:,2)];
        figure(hf3);
        plot3(xyzb(:,1),xyzb(:,2),xyzb(:,3),lt(4,:),'LineWidth',1);
        plot3(xyzc(:,1),xyzc(:,2),xyzc(:,3),lt(3,:),'LineWidth',1);
        if any(rsl1==slk)||any(rsl2==slk)
          if any(rsl1==slk)
            icmp = 1;
          else
            icmp = 2;
          end
          rctrz = mean(xyzb(:,3));
          rctr = [rctrx3(icmp) rslc3(icmp) rctrz];
          dy = xyzb(1,2)-rctr(:,2);
          dy2 = dy*dy;
          dx = sqrt(rrad2-dy2);
          if any(rsl1==slk)
            id1 = id1+1;
            sllim1(id1,:) = [rctrx3(icmp)-dx rctrx3(icmp)+dx];
          else
            id2 = id2+1;
            sllim2(id2,:) = [rctrx3(icmp)-dx rctrx3(icmp)+dx];
          end
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
          if icmp==1&&ifirst1
            plot3(rctr(:,1),rctr(:,2),rctr(:,3),'ms','LineWidth',1);
            [xp,yp,zp] = cylinder(repmat(rrad,4,1),72);
            xp1 = xp+rctr(1);
            yp1 = yp+rctr(2);
            zp1 = -20.0*zp+rctr(3);
%
            if icmprt(k)>1
              plot3(xp1,yp1,zp1,'g-','Color',[0 0.6 0],'LineWidth',1);
              plot3(xp1',yp1',zp1','g-','Color',[0 0.6 0], ...
                    'LineWidth',1);
%
              plot3(xyzm(:,1),xyzm(:,2),xyzm(:,3),'ks');   % ROI center
              plot3(xyzms(:,1),xyzms(:,2),xyzms(:,3),'co');     % Contact points
              xp = xp+xyzm(1);
              yp = yp+xyzm(2);
              zp = -20.0*zp+xyzm(3);
              plot3(xp,yp,zp,'g-','Color',[0 0.8 0],'LineWidth',1);
              plot3(xp',yp',zp','g-','Color',[0 0.8 0],'LineWidth',1);
            else
              plot3(xp1,yp1,zp1,'b-','Color',[0 0 0.75],'LineWidth',1);
              plot3(xp1',yp1',zp1','b-','Color',[0 0 0.75], ...
                    'LineWidth',1);
%
              plot3(xyzl(:,1),xyzl(:,2),xyzl(:,3),'ks');   % ROI center
              plot3(xyzls(:,1),xyzls(:,2),xyzls(:,3),'co');   % Contact points
              xp = xp+xyzl(1);
              yp = yp+xyzl(2);
              zp = -20.0*zp+xyzl(3);
              plot3(xp,yp,zp,'b-','Color',[0 0.5 1],'LineWidth',1);
              plot3(xp',yp',zp','b-','Color',[0 0.5 1],'LineWidth',1);
            end
            ifirst1 = false;
          end
%
          if icmp==2&&ifirst2
            plot3(rctr(:,1),rctr(:,2),rctr(:,3),'ms','LineWidth',1);
            [xp,yp,zp] = cylinder(repmat(rrad,4,1),72);
            xp1 = xp+rctr(1);
            yp1 = yp+rctr(2);
            zp1 = -20.0*zp+rctr(3);
%
            if icmprt(k)>1
              plot3(xp1,yp1,zp1,'g-','Color',[0 0.6 0],'LineWidth',1);
              plot3(xp1',yp1',zp1','g-','Color',[0 0.6 0], ...
                    'LineWidth',1);
%
              plot3(xyzm(:,1),xyzm(:,2),xyzm(:,3),'ks');   % ROI center
              plot3(xyzms(:,1),xyzms(:,2),xyzms(:,3),'co');     % Contact points
              xp = xp+xyzm(1);
              yp = yp+xyzm(2);
              zp = -20.0*zp+xyzm(3);
              plot3(xp,yp,zp,'g-','Color',[0 0.8 0],'LineWidth',1);
              plot3(xp',yp',zp','g-','Color',[0 0.8 0],'LineWidth',1);
            else
              plot3(xp1,yp1,zp1,'b-','Color',[0 0 0.75],'LineWidth',1);
              plot3(xp1',yp1',zp1','b-','Color',[0 0 0.75], ...
                    'LineWidth',1);
%
              plot3(xyzl(:,1),xyzl(:,2),xyzl(:,3),'ks');   % ROI center
              plot3(xyzls(:,1),xyzls(:,2),xyzls(:,3),'co');     % Contact points
              xp = xp+xyzl(1);
              yp = yp+xyzl(2);
              zp = -20.0*zp+xyzl(3);
              plot3(xp,yp,zp,'b-','Color',[0 0.5 1],'LineWidth',1);
              plot3(xp',yp',zp','b-','Color',[0 0.5 1],'LineWidth',1);
            end
            ifirst2 = false;
          end           % End of if region 2 and first slice in region2
        end             % End of if slice in region 1 or 2
%
        if any(rsll==slk)||any(rslm==slk)
          if any(rsll==slk)
            dy = xyzb(1,2)-xyzl(:,2);
            dy2 = dy*dy;
            dx = sqrt(rrad2-dy2);
            idl = idl+1;
            slliml(idl,:) = [xyzl(:,1)-dx xyzl(:,1)+dx];
          else
            dy = xyzb(1,2)-xyzm(:,2);
            dy2 = dy*dy;
            dx = sqrt(rrad2-dy2);
            idm = idm+1;
            sllimm(idm,:) = [xyzm(:,1)-dx xyzm(:,1)+dx];
          end
        end             % End of if slice in lateral or medial region
%
      end               % End of if tibia bone
   end                  % End of slice loop - k
%
   set(gca,'ZDir','reverse');
   view(3);
   axis equal;
   title({dirstr; [fs ' - T1\rho']; 'Blue - Lateral/Green - Medial'}, ...
          'FontSize',16,'FontWeight','bold');
   print('-dpsc2','-r600','-fillpage',pnam3);
   view(-90,90);
   print('-dpsc2','-r600','-fillpage','-append',pnam3);
%
% Get Region 1 Masks
%
   sllim1 = fix(sllim1./pspcs(1));
   maskfr1 = false(npxt,2,nrsl1);      % Mask for region 1 femoral cartilage
   masktr1 = false(npxt,2,nrsl1);      % Mask for region 1 tibial cartilage
%
   for k = 1:nrsl1
%
      slk = rsl1(k);
      maskr1 = false(iszs);
      maskr1(:,sllim1(k,1):sllim1(:,2)) = true;
      maskr1 = repmat(maskr1(:),1,2);
      id1 = slk==rsl;
      maskfr1(:,:,k) = maskf(:,:,id1)&maskr1;
      masktr1(:,:,k) = maskt(:,:,id1)&maskr1;
%
   end
%
% Get Region 2 Masks
%
   sllim2 = fix(sllim2./pspcs(1));
   maskfr2 = false(npxt,2,nrsl2);      % Mask for region 2 femoral cartilage
   masktr2 = false(npxt,2,nrsl2);      % Mask for region 2 tibial cartilage
%
   for k = 1:nrsl2
%
      slk = rsl2(k);
      maskr2 = false(iszs);
      maskr2(:,sllim2(k,1):sllim2(:,2)) = true;
      maskr2 = repmat(maskr2(:),1,2);
      id2 = slk==rsl;
      maskfr2(:,:,k) = maskf(:,:,id2)&maskr2;
      masktr2(:,:,k) = maskt(:,:,id2)&maskr2;
%
   end
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
   rslp = unique([rsl1; rsl2; rsll; rslm]);
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
% Plot 20%/80% ROIs on the Top Half of the Page
%
      if any(rsl1==slk)||any(rsl2==slk)
        if any(rsl1==slk)
          idp = rsl1==slk;
          mask1 = maskfr1(:,1,idp)|masktr1(:,1,idp);  % Superficial cartilage mask
          mask2 = maskfr1(:,2,idp)|masktr1(:,2,idp);  % Deep cartilage mask
          dcmx = 16*cmx/128;
          img1(mask1) = dcmx;  % Blue - Superficial
          img1(mask2) = cmx-dcmx;   % Red - Deep
        else
          idp = rsl2==slk;
          mask1 = maskfr2(:,1,idp)|masktr2(:,1,idp);  % Superficial cartilage mask
          mask2 = maskfr2(:,2,idp)|masktr2(:,2,idp);  % Deep cartilage mask
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
            ' Compartment']; '20%/80% ROIs'}, ...
            'FontSize',12,'FontWeight','bold');
%
% Plot Contact ROIs on the Bottom Half of the Page
%
      img1 = img-cmx-1;
%
      subplot(2,1,2);
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
      title('Contact Points ROIs','FontSize',12,'FontWeight','bold');
%
      print('-dpsc2','-r600','-fillpage','-append',pnam3);
%
   end
%
% Save Masks, ROIS and Slice Information into MAT File
% Note maskp and p are empty/false.
%
   savnam = [mnam(1:end-4) '_rois7.mat'];
   save(savnam,'f','ibone','icmprt','irsl1','irsl2','maskf', ...
               'maskfr1','maskfr2','maskfrl','maskfrm','maskt', ...
               'masktr1','masktr2','masktrl','masktrm','brois', ...
               'rsl','rsl1','rsl2','rsll','rslm','t');
%
   close all;
%
end                     % End of m loop - MAT file loop
%
end                     % End of if skip RHO
%
return
%
% T2star
%
rdir = 'T2S';           % Directory for T2* segmentations
id5 = repmat(3,4,1);    % Default echo time for segmentations
thresh = 1.75;          % Four (4) pixels = 4*0.4375 = 1.75
% threshs = num2str(thresh,'%.3f');      % Threshold as a string
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
for m = 1:nmat
% for m = 4:nmat
% for m = 2:2
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
   else
     leg = 'R';
   end
%
   if contains(st,'Load','IgnoreCase',true)
     ld = 'LD';
   else
     ld = 'UL';
   end
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
%
   rslmn = min(rslb);   % Tibia bone minimum slice
   rslmx = max(rslb);   % Tibia bone maximum slice
%
% Get Slices Within RRAD of the 20% and 80% Width of Tibia Bone
%
   rslctrs = rslmn+(rslmx-rslmn).*[0.2 0.8];     % 1st column = 20%, 2nd column = 80%
   rslc3 = rslctrs*2.0; % Slice thickness = 2 mm
   rsl3 = 2.0*rsl;      % Slice thickness = 2 mm
%
   idx = (rsl3>rslc3-rrad)&(rsl3<rslc3+rrad);    % 1st column = 20%, 2nd column = 80%
   rsl1 = rsl(idx(:,1));               % First compartment slices at 20%
   rsl2 = rsl(idx(:,2));               % Second compartment slices at 80%
   regpx = zeros(2);    % Minimum and maximum pixel coordinates for both lateral/medial regions
   regpx(1,:) = iszs(1);     % Minimums
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_ROIs1.ps'];           % ROI lines print file name
   pnam2 = [fs '_ROIs2.ps'];           % ROI areas print file name
   pnam3 = [fs '_ROIs3.ps'];           % ROI regions print file name
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
% Check for Region Slices
%
      if any(rsl1==slk)||any(rsl2==slk)
        if any(rsl1==slk)
          icol = 1;
        else
          icol = 2;
        end
        dat2 = t{2,k}(:,1);
        xmn = min(dat2);
        xmx = max(dat2);
        if xmn<regpx(1,icol)
          regpx(1,icol) = xmn;
        end
        if xmx>regpx(2,icol)
          regpx(2,icol) = xmx;
        end
      end
%
   end                  % End of k loop - tibia slices
%
% Close Slice Plots
%
   close all;
%
% Get Centers of Regions
%
   rctrx = regpx(1,:)'+diff(regpx)'/2;
   rctrx3 = rctrx*pspcs(1);
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
% Get Compartment for Each Center
%
   ksl = rsl==rsl1(1);
   irsl1 = icmprt(ksl);
   ksl = rsl==rsl2(end);
   irsl2 = icmprt(ksl);
%
% Loop through the Slices to Define Analysis Regions
%
   id1 = 0;
   id2 = 0;
   idl = 0;
   idm = 0;
   ifirst1 = true;
   ifirst2 = true;
%
   nrsl1 = size(rsl1,1);
   nrsl2 = size(rsl2,1);
   nrsll = size(rsll,1);     % Number of lateral slices in contact ROI
   nrslm = size(rslm,1);     % Number of medial slices in contact ROI
%
   sllim1 = zeros(nrsl1,2);  % Column 1 - minimum, column 2 - maximum
   sllim2 = zeros(nrsl1,2);  % Column 1 - minimum, column 2 - maximum
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
        if any(rsl1==slk)||any(rsl2==slk)
          if any(rsl1==slk)
            icmp = 1;
          else
            icmp = 2;
          end
          rctrz = mean(xyzb(:,3));
          rctr = [rctrx3(icmp) rslc3(icmp) rctrz];
          dy = xyzb(1,2)-rctr(:,2);
          dy2 = dy*dy;
          dx = sqrt(rrad2-dy2);
          if any(rsl1==slk)
            id1 = id1+1;
            sllim1(id1,:) = [rctrx3(icmp)-dx rctrx3(icmp)+dx];
          else
            id2 = id2+1;
            sllim2(id2,:) = [rctrx3(icmp)-dx rctrx3(icmp)+dx];
          end
%
          xyzcf = f{1,k}*pspcs(1);     % Assumes square pixels
          npts = size(xyzcf,1);
          if npts>0
            xyzcf = [xyzcf(:,1) repmat(2.0*slk,npts,1) xyzcf(:,2)];
            plot3(xyzcf(:,1),xyzcf(:,2),xyzcf(:,3),lt(1,:), ...
                  'LineWidth',1);
          end
%
          xyzbf = f{2,k}*pspcs(1);     % Assumes square pixels
          npts = size(xyzbf,1);
          xyzbf = [xyzbf(:,1) repmat(2.0*slk,npts,1) xyzbf(:,2)];
          plot3(xyzbf(:,1),xyzbf(:,2),xyzbf(:,3),lt(2,:), ...
                'LineWidth',1);
%
          if icmp==1&&ifirst1
            plot3(rctr(:,1),rctr(:,2),rctr(:,3),'ms','LineWidth',1);
            [xp,yp,zp] = cylinder(repmat(rrad,4,1),72);
            xp1 = xp+rctr(1);
            yp1 = yp+rctr(2);
            zp1 = -18.0*zp+rctr(3);
%
            if icmprt(k)>1
              plot3(xp1,yp1,zp1,'g-','Color',[0 0.6 0],'LineWidth',1);
              plot3(xp1',yp1',zp1','g-','Color',[0 0.6 0], ...
                    'LineWidth',1);
%
              plot3(xyzm(:,1),xyzm(:,2),xyzm(:,3),'ks');   % ROI center
              plot3(xyzms(:,1),xyzms(:,2),xyzms(:,3),'co');     % Contact points
              xp = xp+xyzm(1);
              yp = yp+xyzm(2);
              zp = -18.0*zp+xyzm(3);
              plot3(xp,yp,zp,'g-','Color',[0 0.8 0],'LineWidth',1);
              plot3(xp',yp',zp','g-','Color',[0 0.8 0],'LineWidth',1);
            else
              plot3(xp1,yp1,zp1,'b-','Color',[0 0 0.75],'LineWidth',1);
              plot3(xp1',yp1',zp1','b-','Color',[0 0 0.75], ...
                    'LineWidth',1);
%
              plot3(xyzl(:,1),xyzl(:,2),xyzl(:,3),'ks');   % ROI center
              plot3(xyzls(:,1),xyzls(:,2),xyzls(:,3),'co');   % Contact points
              xp = xp+xyzl(1);
              yp = yp+xyzl(2);
              zp = -18.0*zp+xyzl(3);
              plot3(xp,yp,zp,'b-','Color',[0 0.5 1],'LineWidth',1);
              plot3(xp',yp',zp','b-','Color',[0 0.5 1],'LineWidth',1);
            end
            ifirst1 = false;
          end
%
          if icmp==2&&ifirst2
            plot3(rctr(:,1),rctr(:,2),rctr(:,3),'ms','LineWidth',1);
            [xp,yp,zp] = cylinder(repmat(rrad,4,1),72);
            xp1 = xp+rctr(1);
            yp1 = yp+rctr(2);
            zp1 = -18.0*zp+rctr(3);
%
            if icmprt(k)>1
              plot3(xp1,yp1,zp1,'g-','Color',[0 0.6 0],'LineWidth',1);
              plot3(xp1',yp1',zp1','g-','Color',[0 0.6 0], ...
                    'LineWidth',1);
%
              plot3(xyzm(:,1),xyzm(:,2),xyzm(:,3),'ks');   % ROI center
              plot3(xyzms(:,1),xyzms(:,2),xyzms(:,3),'co');     % Contact points
              xp = xp+xyzm(1);
              yp = yp+xyzm(2);
              zp = -20.0*zp+xyzm(3);
              plot3(xp,yp,zp,'g-','Color',[0 0.8 0],'LineWidth',1);
              plot3(xp',yp',zp','g-','Color',[0 0.8 0],'LineWidth',1);
            else
              plot3(xp1,yp1,zp1,'b-','Color',[0 0 0.75],'LineWidth',1);
              plot3(xp1',yp1',zp1','b-','Color',[0 0 0.75], ...
                    'LineWidth',1);
%
              plot3(xyzl(:,1),xyzl(:,2),xyzl(:,3),'ks');   % ROI center
              plot3(xyzls(:,1),xyzls(:,2),xyzls(:,3),'co');     % Contact points
              xp = xp+xyzl(1);
              yp = yp+xyzl(2);
              zp = -20.0*zp+xyzl(3);
              plot3(xp,yp,zp,'b-','Color',[0 0.5 1],'LineWidth',1);
              plot3(xp',yp',zp','b-','Color',[0 0.5 1],'LineWidth',1);
            end
            ifirst2 = false;
          end           % End of if region 2 and first slice in region2
        end             % End of if slice in region 1 or 2
%
        if any(rsll==slk)||any(rslm==slk)
          if any(rsll==slk)
            dy = xyzb(1,2)-xyzl(:,2);
            dy2 = dy*dy;
            dx = sqrt(rrad2-dy2);
            idl = idl+1;
            slliml(idl,:) = [xyzl(:,1)-dx xyzl(:,1)+dx];
          else
            dy = xyzb(1,2)-xyzm(:,2);
            dy2 = dy*dy;
            dx = sqrt(rrad2-dy2);
            idm = idm+1;
            sllimm(idm,:) = [xyzm(:,1)-dx xyzm(:,1)+dx];
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
% Get Region 1 Masks
%
   sllim1 = fix(sllim1./pspcs(1));
   maskfr1 = false(npxt,2,nrsl1);      % Mask for region 1 femoral cartilage
   masktr1 = false(npxt,2,nrsl1);      % Mask for region 1 tibial cartilage
%
   for k = 1:nrsl1
%
      slk = rsl1(k);
      maskr1 = false(iszs);
      maskr1(:,sllim1(k,1):sllim1(:,2)) = true;
      maskr1 = repmat(maskr1(:),1,2);
      id1 = slk==rsl;
      maskfr1(:,:,k) = maskf(:,:,id1)&maskr1;
      masktr1(:,:,k) = maskt(:,:,id1)&maskr1;
%
   end
%
% Get Region 2 Masks
%
   sllim2 = fix(sllim2./pspcs(1));
   maskfr2 = false(npxt,2,nrsl1);      % Mask for region 2 femoral cartilage
   masktr2 = false(npxt,2,nrsl1);      % Mask for region 2 tibial cartilage
%
   for k = 1:nrsl2
%
      slk = rsl2(k);
      maskr2 = false(iszs);
      maskr2(:,sllim2(k,1):sllim2(:,2)) = true;
      maskr2 = repmat(maskr2(:),1,2);
      id2 = slk==rsl;
      maskfr2(:,:,k) = maskf(:,:,id2)&maskr2;
      masktr2(:,:,k) = maskt(:,:,id2)&maskr2;
%
   end
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
   rslp = unique([rsl1; rsl2; rsll; rslm]);
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
% Plot 20%/80% ROIs on the Top Half of the Page
%
      if any(rsl1==slk)||any(rsl2==slk)
        if any(rsl1==slk)
          idp = rsl1==slk;
          mask1 = maskfr1(:,1,idp)|masktr1(:,1,idp);  % Superficial cartilage mask
          mask2 = maskfr1(:,2,idp)|masktr1(:,2,idp);  % Deep cartilage mask
          dcmx = 16*cmx/128;
          img1(mask1) = dcmx;  % Blue - Superficial
          img1(mask2) = cmx-dcmx;   % Red - Deep
        else
          idp = rsl2==slk;
          mask1 = maskfr2(:,1,idp)|masktr2(:,1,idp);  % Superficial cartilage mask
          mask2 = maskfr2(:,2,idp)|masktr2(:,2,idp);  % Deep cartilage mask
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
            ' Compartment']; '20%/80% ROIs'}, ...
            'FontSize',12,'FontWeight','bold');
%
% Plot Contact ROIs on the Bottom Half of the Page
%
      img1 = img-cmx-1;
%
      subplot(2,1,2);
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
      title('Contact Points ROIs','FontSize',12,'FontWeight','bold');
%
      print('-dpsc2','-r600','-fillpage','-append',pnam3);
%
   end
%
% Save Masks, ROIs and Slice Information into MAT File
% Note maskp and p are empty/false.
%
   savnam = [mnam(1:end-4) '_rois.mat'];
   save(savnam,'f','ibone','icmprt','irsl1','irsl2','maskf', ...
               'maskfr1','maskfr2','maskfrl','maskfrm','maskt', ...
               'masktr1','masktr2','masktrl','masktrm','brois', ...
               'rsl','rsl1','rsl2','rsll','rslm','t');
%
   close all;
%
end                     % End of m loop - MAT file loop
%
return