%#######################################################################
%
%      * ALL SEGmentation as Regions of Interest (ROIs) Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     CSV files to create masks for the cartilage regions of interest
%     in the lateral and medial tibial compartments and lateral and
%     medial femoral compartments.  The masks are saved in MAT files
%     with the series number and ending in "_arois.mat."
%
%     NOTES:  1.  The registered MRI MAT files must be in subject
%             directories starting with "MRIR" and either "Visit1" or
%             "Visit2" visit subdirectories.
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  See rd_dicom.m.
%
%             3.  Femoral trochlear slices without corresponding tibial
%             slices are considered true trochlear slices and are not
%             considered parts of the regions of interest.
%
%             4.  The M-filesv cr_maskf.m, cr_maskt.m, fix_gap.m,
%             in_tri2d.m, lsect2.m, lsect2a.m, mk_tri_2d.m,
%             mk_tri_2df.m, near2.m, rd_rois3.m, and rd_roi6.m must be
%             in the current directory or path.
%
%     04-Nov-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Setup Color Maps, Line Types and Color, and Legend Strings
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
% rmap = [[0.85 0 0]; [0 0 0.85]; [0 0.7 0]; gmap];     % Red/blue/green color map for ROIs
%
lt = ['g.-'; 'b.-'; 'c.-'; 'r.-'; 'y.-'; 'm.-']; % Line color and type
legds = ['FEM_CART'; 'FEM_BONE'; 'TIB_CART'; 'TIB_BONE'];
%
% Segmentation Parameter and Trochlear Region
%
dist = 7.5;             % Half distance from cartilage to bone endpoints
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
% Get Subject and Visit Indices
%
subj = eval(dirstr(1:2));
vnum = eval(dirstr(end));
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
  error(' *** ERROR in all_seg_rois:  Not four (4) MAT files!')
end
%
% Loop through T1rho MAT Files
%
for m = 1:nmat
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
% Read Segmentations
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
   rslt = intersect(rslc,rslb);        % Ensure unique slices in sorted order
   nrslt = size(rslt,1);
%
% Get Femoral Trochlear Slices with Corresponding Tibial Slices
%
   rslftc = brois(1).rois(1).roi(3).imageno';    % Trochlear cartilage
   rslftb = brois(1).rois(2).roi(3).imageno';    % Trochlear bone
   rslft = intersect(rslftc,rslftb);   % Ensure unique slices in sorted order
   sl2rm = setdiff(rslft,rslt);        % Slices to remove (true trochlea)
   rslf = setdiff(rslf,sl2rm);         % Remove trochlear slices
   nrslf = size(rslf,1);
%
% Both Femur and Tibia Slices
%
   rsl = union(rslf,rslt);
   nrsl = size(rsl,1);
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_arois1.ps'];          % ROI lines print file name
   pnam2 = [fs '_arois2.ps'];          % ROI areas print file name
%
   f = cell(2,nrsl);    % Femur coordinates (1 - cartilage, 2 - bone)
   t = cell(2,nrsl);    % Tibia coordinates (1 - cartilage, 2 - bone)
%
   ibone = false(nrsl,2);              % Bone exists? (columns:  1 - femur, 2 - tibia)
   icmprtf = zeros(nrsl,1);            % Femoral compartments (1 - lateral, 2 - medial)
   icmprtt = zeros(nrsl,1);            % Tibial compartments (1 - lateral, 2 - medial)
   npxt = prod(iszs);                  % Total number of pixels in a slice
   maskf = false(npxt,nrsl);           % Mask for femoral cartilage
   maskt = false(npxt,nrsl);           % Mask for tibial cartilage
%
% Loop through All Slices
%
   for k = 1:nrsl       % Loop through all slices
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
                   dat = fix_gap(dat);
                   if lb==1
                     f{lc,k} = dat;    % Femur
                     icmprtf(k) = n;
                   else
                     t{lc,k} = dat;    % Tibia
                     icmprtt(k) = n;
                   end
                   lh(idxl) = plot(dat(:,1),dat(:,2),lt(idxl,:));
                   idl(idxl) = true;
                 end
              end       % End of n loop - lateral/medial/trochlear loop
            end         % End of if - any rois for this slice
         end            % End of lc loop - cartilage/bone loop
         if lb==1
           ibone(k,lb) = all(~cellfun('isempty',f(:,k)));
         else
           ibone(k,lb) = all(~cellfun('isempty',t(:,k)));
         end
      end               % End of lb loop - femur/tibia loop
%
% Assign Trochlear Slices to Medial or Lateral Compartment
%
      if icmprtf(k)==3;
        icmprtf(k) = icmprtt(k);       % Use tibia to assign compartment
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
        maskf(:,k) = cr_maskf(f(:,k),iszs,dist);
      end
%
      if ibone(k,2)
        maskt(:,k) = cr_maskt(t(:,k),iszs);
      end
%
% Plot ROIs
%
      mask = maskf(:,k)|maskt(:,k);    % All cartilage mask
      cmx = max(img(:));
      img1 = img-cmx-1;
      idmsk = ~mask;
      img2 = img;
      img2(idmsk) = img(idmsk)-max(img(:))-1;
      dcmx = 16*cmx/128;
      img1(maskf(:,k)) = dcmx;         % Blue - Femur
      img1(maskt(:,k)) = cmx-dcmx;     % Red - Tibia
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
% Save Masks, ROIs and Slice Information into MAT File
%
   savnam = [mnam(1:end-4) '_arois.mat'];
   save(savnam,'brois','f','ibone','icmprtf','icmprtt','maskf', ...
               'maskt','nrsl','nrslf','nrslt','rsl','rslf','rslt','t');
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
  error(' *** ERROR in all_seg_rois:  Not four (4) MAT files!')
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
   rslt = intersect(rslc,rslb);        % Ensure unique slices in sorted order
   nrslt = size(rslt,1);
%
% Get Femoral Trochlear Slices with Corresponding Tibial Slices
%
   rslftc = brois(1).rois(1).roi(3).imageno';    % Trochlear cartilage
   rslftb = brois(1).rois(2).roi(3).imageno';    % Trochlear bone
   rslft = intersect(rslftc,rslftb);   % Ensure unique slices in sorted order
   sl2rm = setdiff(rslft,rslt);        % Slices to remove (true trochlea)
   rslf = setdiff(rslf,sl2rm);         % Remove trochlear slices
   nrslf = size(rslf,1);
%
% Both Femur and Tibia Slices
%
   rsl = union(rslf,rslt);
   nrsl = size(rsl,1);
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_arois1.ps'];          % ROI lines print file name
   pnam2 = [fs '_arois2.ps'];          % ROI areas print file name
%
   f = cell(2,nrsl);    % Femur coordinates (1 - cartilage, 2 - bone)
   t = cell(2,nrsl);    % Tibia coordinates (1 - cartilage, 2 - bone)
%
   ibone = false(nrsl,2);              % Bone exists? (1 - femur, 2 - tibia)
   icmprtf = zeros(nrsl,1);            % Femoral compartments (1 - lateral, 2 - medial)
   icmprtt = zeros(nrsl,1);            % Tibal compartments (1 - lateral, 2 - medial)
   npxt = prod(iszs);                  % Total number of pixels in a slice
   maskf = false(npxt,2,nrsl);         % Mask for femoral cartilage
   maskt = false(npxt,2,nrsl);         % Mask for tibial cartilage
%
% Loop through All Slices
%
   for k = 1:nrsl       % Loop through all slices
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
              for n = 1:nc   % Loop through compartments - lateral = 1, medial = 2 and trochlea = 3
                 idxs = brois(lb).rois(lc).roi(n).imageno==slk;
                 if any(idxs)
                   dat = cell2mat( ...
                                 brois(lb).rois(lc).roi(n).data(idxs)');
                   dat = fix_gap(dat);
                   if lb==1
                     f{lc,k} = dat;    % Femur
                     icmprtf(k) = n;
                   else
                     t{lc,k} = dat;    % Tibia
                     icmprtt(k) = n;
                   end
                   lh(idxl) = plot(dat(:,1),dat(:,2),lt(idxl,:));
                   idl(idxl) = true;
                 end
              end       % End of n loop - lateral/medial/trochlear loop
            end         % End of if - any rois for this slice
         end            % End of lc loop - cartilage/bone loop
         if lb==1
           ibone(k,lb) = all(~cellfun('isempty',f(:,k)));
         else
           ibone(k,lb) = all(~cellfun('isempty',t(:,k)));
         end
      end               % End of lb loop - femur/tibia loop
%
% Assign Trochlear Slices to Medial or Lateral Compartment
%
      if icmprtf(k)==3;
        icmprtf(k) = icmprtt(k);       % Use tibia to assign compartment
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
        maskf(:,k) = cr_maskf(f(:,k),iszs,dist);
      end
%
      if ibone(k,2)
        maskt(:,k) = cr_maskt(t(:,k),iszs);
      end
%
% Plot ROIs
%
      mask = maskf(:,k)|maskt(:,k);    % All cartilage mask
      cmx = max(img(:));
      img1 = img-cmx-1;
      idmsk = ~mask;
      img2 = img;
      img2(idmsk) = img(idmsk)-max(img(:))-1;
      dcmx = 16*cmx/128;
      img1(maskf(:,k)) = dcmx;         % Blue - Femur
      img1(maskt(:,k)) = cmx-dcmx;     % Red - Tibia
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
% Save Masks, ROIs and Slice Information into MAT File
%
   savnam = [mnam(1:end-4) '_arois.mat'];
   save(savnam,'brois','f','ibone','icmprtf','icmprtt','maskf', ...
               'maskt','nrsl','nrslf','nrslt','rsl','rslf','rslt','t');
%
   close all;
%
end                     % End of m loop - MAT file loop
%
return