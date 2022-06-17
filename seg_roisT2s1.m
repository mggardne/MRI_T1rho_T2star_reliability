%#######################################################################
%
%     * SEGmentation to Regions of Interest (ROIs) T2* 1 Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     CSV files to create masks for cylindrical regions of interest in
%     the lateral and medial tibial compartments.  The masks are saved
%     in MAT files with the series number and ending in "_rois.mat."
%
%     Segmentations overlaid on the first echo time.
%
%     NOTES:  1.  For use with subject 08-AS on Visit 2 for the T2* for
%             the unloaded left leg.
%
%             2.  The registered MRI MAT files must be in subject
%             directories starting with "MRIR" and either "Visit1" or
%             "Visit2" visit subdirectories.
%
%             3.  T2* MAT files must start with "T2star_S".  See
%             rd_dicom.m or rd_dicomT2s.
%
%             4.  M-file circ_plt.m, rd_rois.m and rd_roi6.m must be in
%             the current directory or path.
%
%     07-Mar-2022 * Mack Gardner-Morse
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
% Directory with MAT Files as a String
%
dirstr = split(pwd,filesep);
vstr = dirstr{end};
vstr = [vstr(1:end-1) ' ' vstr(end)];
dirstr = dirstr{end-1};
dirstr = split(dirstr,' ');
dirstr = [dirstr{end} ' - ' vstr];
%
% T2star
%
rdir = 'T2S';           % Directory for T2* segmentations
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
for m = nmat:nmat
% for m = 2:2
%
   mnam = mnams{m};
   load(mnam);
   fs = ['S' snt];      % Series number prefaced with a 'S'
   id5m = 1;            % Use first echo time
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
% Get Tibia Bone Slices
%
   rsl = brois(2).rois(2).slice;
   rsl = unique(rsl);        % Ensure unique slices in sorted order
   nrsl = size(rsl,1);
   rslmn = rsl(1);
   rslmx = rsl(nrsl);
%
% Get Slices Within RRAD of the 20% and 80% Width of Tibia Bone
%
   rslctrs = rslmn+(rslmx-rslmn).*[0.2 0.8];     % 1st column = 20%, 2nd column = 80%
   rslc3 = rslctrs*2.0; % Slice thickness = 2 mm
   rsl3 = 2.0*rsl;      % Slice thickness = 2 mm
%
   idx = (rsl3>rslc3-rrad)&(rsl3<rslc3+rrad);    % 1st column = 20%, 2nd column = 80%
   rsl1 = rsl(idx(:,1));
   rsl2 = rsl(idx(:,2));
   regpx = zeros(2);    % Minimum and maximum pixel coordinates for both lateral/medial regions
   regpx(1,:) = iszs(1);     % Minimums
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_ROIs1_1.ps'];         % ROI lines print file name
   pnam2 = [fs '_ROIs2_1.ps'];         % ROI areas print file name
   pnam3 = [fs '_ROIs3_1.ps'];         % ROI regions print file name
%
   f = cell(2,nrsl);    % Femur coordinates (1 - cartilage, 2 - bone)
   t = cell(2,nrsl);    % Tibia coordinates (1 - cartilage, 2 - bone)
%
   ibone = false(nrsl,2);              % false - femur, true - tibia
   icmprt = ones(nrsl,1);              % 1 - lateral, 2 - medial, 3 - trochlea
   npxt = prod(iszs);                  % Total number of pixels in a slice
   maskf = false(npxt,2,nrsl);         % Mask for femoral cartilage
   maskt = false(npxt,2,nrsl);         % Mask for tibial cartilage
%
% Loop through Tibia Bones Slices
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
                   icmprt(k) = n;
                   dat = cell2mat(brois(lb).rois(lc).roi(n).data(idxs)');
                   if lb==1
                     f{lc,k} = dat;  % Femur
                   else
                     t{lc,k} = dat;  % Tibia
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
      if any(rsl1==slk|rsl2==slk)
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
   end
%
% Get Centers of Regions
%
   rctrx = regpx(1,:)'+diff(regpx)'/2;
   rctrx3 = rctrx*pspcs(1);
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
   ifirst1 = true;
   ifirst2 = true;
   nrsl1 = size(rsl1,1);
   nrsl2 = size(rsl2,1);
   sllim1 = zeros(nrsl1,2);  % Column 1 - minimum, column 2 - maximum
   sllim2 = zeros(nrsl1,2);  % Column 1 - minimum, column 2 - maximum
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
        if any(rsl1==slk|rsl2==slk)
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
          plot3(rctr(:,1),rctr(:,2),rctr(:,3),'ms','LineWidth',1);
          [xp,yp] = circ_plt(rrad,rctr([1 2]));
          if icmp==1&&ifirst1
%            [xp,yp,zp] = cylinder(repmat(rrad,4,1),72);
            [xp,yp,zp] = cylinder(repmat(rrad,2,1),72);
            xp = xp+rctr(1);
            yp = yp+rctr(2);
            zp = -15.0*zp+rctr(3)+2;
%            zp = -12.0*zp+rctr(3)+2;

%            plot3(xp,yp,zp,'k-','LineWidth',1);
%            plot3(xp',yp',zp','k-','LineWidth',1);

            if icmprt(k)>1
              plot3(xp,yp,zp,'g-','Color',[0 0.7 0],'LineWidth',1);
              plot3(xp',yp',zp','g-','Color',[0 0.7 0],'LineWidth',1);
            else
              plot3(xp,yp,zp,'b-','Color',[0 0 0.8],'LineWidth',1);
              plot3(xp',yp',zp','b-','Color',[0 0 0.8],'LineWidth',1);
            end
            ifirst1 = false;
          end
          if icmp==2&&ifirst2
%            [xp,yp,zp] = cylinder(repmat(rrad,4,1),72);
            [xp,yp,zp] = cylinder(repmat(rrad,2,1),72);
            xp = xp+rctr(1);
            yp = yp+rctr(2);
            zp = -15.0*zp+rctr(3)+2;
%            zp = -12.0*zp+rctr(3)+2;

%            plot3(xp,yp,zp,'k-','LineWidth',1);
%            plot3(xp',yp',zp','k-','LineWidth',1);

            if icmprt(k)>1
              plot3(xp,yp,zp,'g-','Color',[0 0.7 0],'LineWidth',1);
              plot3(xp',yp',zp','g-','Color',[0 0.7 0],'LineWidth',1);
            else
              plot3(xp,yp,zp,'b-','Color',[0 0 0.8],'LineWidth',1);
              plot3(xp',yp',zp','b-','Color',[0 0 0.8],'LineWidth',1);
            end
            ifirst2 = false;
          end
        end
      end
   end
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
      slk = rsl1(k);
      maskr1 = false(iszs);
      maskr1(:,sllim1(k,1):sllim1(:,2)) = true;
      maskr1 = repmat(maskr1(:),1,2);
      id1 = slk==rsl;
      maskfr1(:,:,k) = maskf(:,:,id1)&maskr1;
      masktr1(:,:,k) = maskt(:,:,id1)&maskr1;
%
      img = squeeze(v(:,:,slk,id5m));
      mask1 = maskfr1(:,1,k)|masktr1(:,1,k);     % Superficial cartilage mask
      mask2 = maskfr1(:,2,k)|masktr1(:,2,k);     % Deep cartilage mask
      cmx = max(img(:));
      img1 = img-cmx-1;
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
      title({[fs ' Slice ' int2str(slk)]; 'Region 1'},'FontSize',16, ...
            'FontWeight','bold');
%
      print('-dpsc2','-r600','-fillpage','-append',pnam3);
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
      img = squeeze(v(:,:,slk,id5m));
      mask1 = maskfr2(:,1,k)|masktr2(:,1,k);     % Superficial cartilage mask
      mask2 = maskfr2(:,2,k)|masktr2(:,2,k);     % Deep cartilage mask
      cmx = max(img(:));
      img1 = img-cmx-1;
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
      title({[fs ' Slice ' int2str(slk)]; 'Region 2'},'FontSize',16, ...
            'FontWeight','bold');
%
      print('-dpsc2','-r600','-fillpage','-append',pnam3);
%
   end
%
% Save Masks, ROIs and Slice Information into MAT File
% Note maskp and p are empty/false.
%
%   savnam = [mnam(1:end-4) '_rois.mat'];
%   save(savnam,'f','ibone','icmprt','irsl1','irsl2','maskf', ...
%               'maskfr1','maskfr2','maskt','masktr1','masktr2', ...
%               'brois','rsl','rsl1','rsl2','t');
%
%   close all;
%
end                     % End of m loop - MAT file loop
%
return