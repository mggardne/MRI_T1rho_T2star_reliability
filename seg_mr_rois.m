%#######################################################################
%
%               * SEGmentation to Meniscus Reliability *
%                * Regions of Interest (ROIs) Program *
%
%          M-File which reads the registered MRI data and meniscus
%     segmentation CSV files to create masks for the meniscus regions
%     of interest in the lateral and medial compartments.  The masks
%     are saved in MAT files with the series number and ending in
%     "_mrois.mat."
%
%     NOTES:  1.  The registered MRI MAT files must be in subdirectories
%             "Visit1" or "Visit2" in subject directories "MRIR *"
%             where "*" is the subject number.
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  See rd_dicom.m.
%
%             3.  M-files cr_maskm.m, fnd_bbox.m, in_tri2d.m,
%             mk_tri_m.m, and rd_roi6hp.m must be in the current
%             directory or path.
%
%             4.  Currently only does T1rho segmentations.
%
%     24-Jun-2024 * Mack Gardner-Morse * First version for MRIR study
%                                        data.
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
lt = ['g.-'; 'b.-'; 'r.-'; 'c.-'; 'y.-'; 'm.-']; % Line color and type
legds = {'ANT_LAT'; 'ANT_MED'; 'POST_LAT'; 'POST_MED'};
nlegds = size(legds,1);
%
ipad = 25;              % Bounding box padding in pixels
%
% Compartment Labels
%
cmprt = {'Lateral'; 'Medial'};
%
% Knee Labels
%
legs = {' Left '; ' Right '};
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
% Loop through Subjects
%
for ks = 1:nsubj
%
% T1rho Series and Segmentations
%
   id5m = 1;            % Use first spin lock time for plots
   rdir = 'RHO';        % Directory for T1rho segmentations
%
% Get Subject Directory, Name and Number
%
   sdir = sdirs{ks};    % Current subject directory
   subj = eval(sdir(6:7));             % Subject number
   if subj<10
     snum = ['00' int2str(subj)];
   else
     snum = ['0' int2str(subj)];
   end
%
% Loop through Visits
%
   for kv = 1:nvdir
%
% Get Visit Directory 
%
      vdir = vdirs{kv};
      svdir = fullfile(sdir,vdir);
%
% Get T1rho MAT Files in Directory
%
      d = dir(fullfile(svdir,'T1rho_S*.mat'));
      mnams = {d.name}';
      idr = contains(mnams,'roi','IgnoreCase',true);  % ROI files
      mnams = mnams(~idr);
      idc = contains(mnams,'chk','IgnoreCase',true);  % Check files
      mnams = mnams(~idc);
      id3 = contains(mnams,'_3_');     % _3_?.mat files
      mnams = mnams(~id3);
      nmat = size(mnams,1);
%
      if nmat~=4
        error(' *** ERROR in seg_mr_rois:  Not four (4) MAT files!')
      end
%
% Loop through T1rho MAT Files
%
      for m = 1:nmat
%
         mnam = mnams{m};
         load(fullfile(svdir,mnam),'iszs','snt','st','v');
         fs = ['S' snt];     % Series number prefaced with a 'S'
%
% Parse Series Text for Leg and Load
%
         leg = upper(st(1));
         if strcmpi(leg,'L')
           lstr = legs{1};
         else
           lstr = legs{2};
         end
%
         if contains(st,'Load','IgnoreCase',true)
           ld = 'LD';   % Loaded
           ldstr = 'Loaded';
         else
           ld = 'UL';   % Unloaded
           ldstr = 'Unloaded';
         end
%
% Get Segmentation File Name and Read ROIs
%
         segdir = fullfile(svdir,rdir,'Meniscus');    % Segmentation path
         rnam = dir(fullfile(segdir,[snum '_' leg ...
                    '_SAG_MENISCI_RHO_' ld '*.csv']));
         rnam = {rnam.name}';
         idx = ~contains(rnam,'_REL','IgnoreCase',true);   % Reliability CSV file
         rnam = rnam(idx);
         idx = contains(rnam,'_EDIT'); % Check for edited files
         if any(idx)
           rnam = rnam(idx);
         end
%
         if size(rnam,1)>1
           error([' *** ERROR in seg_mr_rois:  More than one T1rho', ...
                  ' segmentation file for subject ' snum ' on the ', ...
                  lstr ' leg in ' ldstr ' condition!']);
         end
%
         rnam = char(rnam);
         rois = rd_roi6hp(fullfile(segdir,rnam),true);
         nrois = size(rois,1);
%
% Get MRI Slices
%
         rsl = [rois.imageno]';
         rsl = unique(rsl);  % Ensure unique slices in sorted order
         nrsl = size(rsl,1);
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
         pnam1 = fullfile(svdir,[fs '_ROIs1.ps']);     % ROI lines print file name
         pnam2 = fullfile(svdir,[fs '_ROIs2.ps']);     % ROI areas print file name
%
         mcoords = cell(nrsl,nrois);   % Meniscus coordinates
%
         iap = zeros(nrsl,nrois);           % -1 - anterior, 1 - posterior
         icmprt = zeros(nrsl,nrois);        % -1 - lateral, 1 - medial
         npxt = prod(iszs);                 % Total number of pixels in a slice
         maskm = false(npxt,nrsl,nrois);    % Masks for all meniscus segmentations
         masks = false(npxt,nrsl);          % Masks for meniscus slices
%
% Loop through Slices
%
         for k = 1:nrsl      % Loop through slices
%
% Get T1 Slice Image
%
            slk = rsl(k);    % Slice number
            img = squeeze(v(:,:,slk,id5m));
%
            figure;
            orient landscape;
            imagesc(img);
            colormap gray;
            axis image;
            axis off;
            hold on;
%
            lh = gobjects(nrois,1);    % Line graphic handles
            idl = false(nrois,1);
            idxs = false(nlegds,1);
%
% Get ROI Data for this Slice and Plot ROIs
%
            for l = 1:nrois            % Loop through ROIs
               idxr = rois(l).imageno==slk;
               if any(idxr)
                 mcoords(k,l) = rois(l).data(idxr);   % Meniscus coordinates
%
                 snam = rois(l).name;  % Get ROI name
%
                 if strcmpi(snam(1),'A')
                   iap(k,l) = -1;      % -1 - anterior
                 else
                   iap(k,l) = 1;       % 1 - posterior
                 end
                 if strcmpi(snam(end-1),'L')
                   icmprt(k,l) = -1;   % -1 - lateral
                 else
                   icmprt(k,l) = 1;    % 1 - medial
                 end
%
% Create Logical Masks for this Slice
%
                 maskm(:,k,l) = cr_maskm(mcoords{k,l},iszs);    % Masks for all segmentations
                 masks(:,k) = masks(:,k)|maskm(:,k,l);          % Masks for slices
%
                 idx = iap(k,l)+icmprt(k,l)/2+2.5;    % Index to legend entry
                 idxs(idx) = true;
%
                 lh(l) = plot(mcoords{k,l}(:,1),mcoords{k,l}(:,2), ...
                              lt(idx,:));
                 idl(l) = true;
               end      % End of if - any rois for this slice
            end         % End of rois loop
%
% Add Legend and Title, and Print Slice Plots
%
            bbx = fnd_bbox(masks(:,k),iszs,ipad);
            axis(bbx);
            legend(lh(idl),legds(idxs,:),'Interpreter','none');
            title({[snum ' Visit ' int2str(kv) ' ' lstr ...
                  'Leg T1\rho'];[ldstr ' ' fs ' Slice ', ...
                  int2str(slk)]; [cmprt{any(icmprt(k,:)==1)+1}, ...
                  ' Compartment']},'FontSize',16,'FontWeight','bold');
            if k==1
              print('-dpsc2','-r600','-fillpage',pnam1);
            else
              print('-dpsc2','-r600','-fillpage','-append',pnam1);
            end
%
% Plot ROIs
%
            cmx = max(img(:));
            img1 = img-cmx-1;
            dcmx = 16*cmx/128;
            maskms = squeeze(maskm(:,k,idl));
            img1(maskms(:,iap(k,idl)==-1)) = dcmx;    % Blue - Anterior
            img1(maskms(:,iap(k,idl)==1)) = cmx-dcmx; % Red - Posterior
%
            figure;
            orient landscape;
            imagesc(img1);
            colormap(cmap);
            caxis([-cmx cmx]);
            axis image;
            axis(bbx);
            axis off;
            title({[snum ' Visit ' int2str(kv) ' ' lstr ...
                  'Leg T1\rho'];[ldstr ' ' fs ' Slice ', ...
                  int2str(slk)]; [cmprt{any(icmprt(k,:)==1)+1}, ...
                  ' Compartment']; ['Blue - Anterior and Red - ', ...
                  'Posterior']},'FontSize',16,'FontWeight','bold');
%
            if k==1
              print('-dpsc2','-r600','-fillpage',pnam2);
            else
              print('-dpsc2','-r600','-fillpage','-append',pnam2);
            end
%
         end            % End of slices loop
%
% Get Bounding Box for All Slices
%
         bbox = fnd_bbox(masks,iszs,ipad);
%
% Get Compartment Identifier for Slices
%
         idl = sum(icmprt,2)<0;        % Index to lateral slices
         rsll = rsl(idl);              % Lateral compartment slices
         nl = size(rsll,1);            % Number of lateral slices
         rslm = rsl(~idl);             % Medial compartment slices
         nm = size(rslm,1);            % Number of medial slices
%
% Get Compartment and AP Specific Masks
%
         maskml = maskm(:,idl,:);      % Masks for lateral slices
         maskmm = maskm(:,~idl,:);     % Masks for medial slices
%
         idapl = sum(iap(idl,:));
         idal = idapl<0;               % Column index to anterior lateral masks
         if any(idal)
           maskmal = maskml(:,:,idal); % Anterior lateral mask
         else
           maskmal = false(npxt,nl);   % Anterior lateral mask
         end
         idpl = idapl>0;               % Column index to posterior lateral masks
         if any(idpl)
           maskmpl = maskml(:,:,idpl); % Posterior lateral mask
         else
           maskmpl = false(npxt,nl);   % Posterior lateral mask
         end
%
         idapm = sum(iap(~idl,:));
         idam = idapm<0;               % Column index to anterior medial masks
         if any(idam)
           maskmam = maskmm(:,:,idam); % Anterior medial mask
         else
           maskmam = false(npxt,nm);   % Anterior medial mask
         end
         idpm = idapm>0;               % Column index to posterior medial masks
         if any(idpm)
           maskmpm = maskmm(:,:,idpm); % Posterior medial mask
         else
           maskmpm = false(npxt,nm);   % Posterior medial mask
         end
%
% Get Compartment and AP Specific Meniscus Coordinates
%
         mcoordl = mcoords(idl,:);     % Meniscus coordinates for lateral slices
         mcoordm = mcoords(~idl,:);    % Meniscus coordinates for medial slices
%
         mcoordal = mcoordl(:,idal);   % Anterior lateral meniscus coordinates
         mcoordpl = mcoordl(:,idpl);   % Posterior lateral meniscus coordinates
%
         mcoordam = mcoordm(:,idam);   % Anterior medial meniscus coordinates
         mcoordpm = mcoordm(:,idpm);   % Posterior medial meniscus coordinates
%
% Save Masks, ROIS and Slice Information into MAT File
%
         savnam = [mnam(1:end-4) '_mrois.mat'];
         savnam = fullfile(svdir,savnam);
         save(savnam,'bbox','cmprt','kv','ld','legs','maskmal', ...
              'maskmam','maskmpl','maskmpm','masks','mcoordal', ...
              'mcoordam','mcoordpl','mcoordpm','rsl','rsll','rslm');
%
         close all;
%
      end               % End of m (MAT file) loop
   end                  % End of kv (visit) loop
end                     % End of ks (subject) loop
%
return
%
% T2star
%
   rdir = 'T2';         % Directory for T2* segmentations
   id5m = 4;            % Default echo time (5 ms) for segmentations
%
% Get T2* MAT Files in Directory
%
   d = dir(fullfile(svdir,'T2star_S*.mat'));
   mnams = {d.name}';
   idr = contains(mnams,'roi','IgnoreCase',true);     % ROI files
   mnams = mnams(~idr);
   idr = contains(mnams,'chk','IgnoreCase',true);     % Check files
   mnams = mnams(~idr);
   nmat = size(mnams,1);
%
% Loop through T2* MAT Files
%
   for m = 1:nmat
%
      mnam = mnams{m};
      load(fullfile(svdir,mnam),'iszs','snt','st','v');
      fs = ['S' snt];   % Series number prefaced with a 'S'
%
% Parse Series Text for Leg
%
      leg = upper(st(1));
      if strcmpi(leg,'L')
        lstr = legs{1};
      else
        lstr = legs{2};
      end
%
% Get Segmentation File Name and Read ROIs
%
      segdir = fullfile(sdir,leg,rdir);     % Segmentation path
      rnam = dir(fullfile(segdir,[sdir '_' leg ...
                 '_sag_menisci_T2S_*.csv']));
      rnam = {rnam.name}';
      idx = ~contains(rnam,'_REL','IgnoreCase',true);
      rnam = rnam(idx);
      idx = contains(rnam,'_EDIT');    % Check for edited files
      if any(idx)
        rnam = rnam(idx);
      end
%
      if size(rnam,1)>1
        error([' *** ERROR in seg_m_rois:  More than one T2*', ...
               ' segmentation file for subject ' sdir ' on the ', ...
               lstr ' leg!']);
      end
%
      rnam = char(rnam);
      rois = rd_roi6(fullfile(segdir,rnam),true);
      nrois = size(rois,1);
%
% Get MRI Slices
%
      rsl = [rois.imageno]';
      rsl = unique(rsl);     % Ensure unique slices in sorted order
      nrsl = size(rsl,1);
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
      pnam1 = fullfile(sdir,[fs '_ROIs1.ps']);   % ROI lines print file name
      pnam2 = fullfile(sdir,[fs '_ROIs2.ps']);   % ROI areas print file name
%
      mcoords = cell(nrsl,nrois);      % Meniscus coordinates
%
      iap = zeros(nrsl,nrois);         % -1 - anterior, 1 - posterior
      icmprt = zeros(nrsl,nrois);      % -1 - lateral, 1 - medial
      npxt = prod(iszs);               % Total number of pixels in a slice
      maskm = false(npxt,nrsl,nrois);  % Masks for all meniscus segmentations
      masks = false(npxt,nrsl);        % Masks for meniscus slices
%
% Loop through Slices
%
      for k = 1:nrsl    % Loop through slices
%
% Get T2 Slice Image
%
         slk = rsl(k);  % Slice number
         img = squeeze(v(:,:,slk,id5m));
%
         figure;
         orient landscape;
         imagesc(img);
         colormap gray;
         axis image;
         axis off;
         hold on;
%
         lh = gobjects(nrois,1);       % Line graphic handles
         idl = false(nrois,1);
         idxs = false(nlegds,1);
%
% Get ROI Data for this Slice and Plot ROIs
%
         for l = 1:nrois               % Loop through ROIs
            idxr = rois(l).imageno==slk;
            if any(idxr)
              mcoords(k,l) = rois(l).data(idxr); % Meniscus coordinates
%        
              snam = rois(l).name;     % Get ROI name
%        
              if strcmpi(snam(1),'A')
                iap(k,l) = -1;         % -1 - anterior
              else
                iap(k,l) = 1;          % 1 - posterior
              end
              if strcmpi(snam(end-1),'L')
                icmprt(k,l) = -1;      % -1 - lateral
              else
                icmprt(k,l) = 1;       % 1 - medial
              end
%
% Create Logical Masks for this Slice
%
              maskm(:,k,l) = cr_maskm(mcoords{k,l},iszs);  % Masks for all segmentations
              masks(:,k) = masks(:,k)|maskm(:,k,l);        % Masks for slices
%
              idx = iap(k,l)+icmprt(k,l)/2+2.5;
              idxs(idx) = true;
%
              lh(l) = plot(mcoords{k,l}(:,1),mcoords{k,l}(:,2), ...
                           lt(idx,:));
              idl(l) = true;
            end         % End of if - any rois for this slice
         end            % End of rois loop
%
% Add Legend and Title, and Print Slice Plots
%
         bbx = fnd_bbox(masks(:,k),iszs,ipad);
         axis(bbx);
         legend(lh(idl),legds(idxs,:),'Interpreter','none');
         title({[sdir lstr 'Leg T2*']; [fs ' Slice ' int2str(slk)]; ...
               [cmprt{any(icmprt(k,:)==1)+1} ' Compartment']}, ...
               'FontSize',16,'FontWeight','bold');
         if k==1
           print('-dpsc2','-r600','-fillpage',pnam1);
         else
           print('-dpsc2','-r600','-fillpage','-append',pnam1);
         end
%
% Plot ROIs
%
         cmx = max(img(:));
         img1 = img-cmx-1;
         dcmx = 16*cmx/128;
         maskms = squeeze(maskm(:,k,idl));
         img1(maskms(:,iap(k,idl)==-1)) = dcmx;       % Blue - Anterior
         img1(maskms(:,iap(k,idl)==1)) = cmx-dcmx;    % Red - Posterior
%
         figure;
         orient landscape;
         imagesc(img1);
         colormap(cmap);
         caxis([-cmx cmx]);
         axis image;
         axis(bbx);
         axis off;
         title({[sdir lstr 'Leg T2*']; [fs ' Slice ' int2str(slk)]; ...
               [cmprt{any(icmprt(k,:)==1)+1} ' Compartment']; ...
               'Blue - Anterior and Red - Posterior'}, ...
               'FontSize',16,'FontWeight','bold');
%
         if k==1
           print('-dpsc2','-r600','-fillpage',pnam2);
         else
           print('-dpsc2','-r600','-fillpage','-append',pnam2);
         end
%
      end               % End of slices loop
%
% Get Bounding Box for All Slices
%
      bbox = fnd_bbox(masks,iszs,ipad);
%
% Get Compartment Identifier for Slices
%
      idl = sum(icmprt,2)<0;           % Index to lateral slices
      rsll = rsl(idl);                 % Lateral compartment slices
      nl = size(rsll,1);               % Number of lateral slices
      rslm = rsl(~idl);                % Medial compartment slices
      nm = size(rslm,1);               % Number of medial slices
%
% Get Compartment and AP Specific Masks
%
      maskml = maskm(:,idl,:);         % Masks for lateral slices
      maskmm = maskm(:,~idl,:);        % Masks for medial slices
%
      idapl = sum(iap(idl,:));
      idal = idapl<0;                  % Column index to anterior lateral masks
      if any(idal)
        maskmal = maskml(:,:,idal);    % Anterior lateral mask
      else
        maskmal = false(npxt,nl);      % Anterior lateral mask
      end
      idpl = idapl>0;                  % Column index to posterior lateral masks
      if any(idpl)
        maskmpl = maskml(:,:,idpl);    % Posterior lateral mask
      else
        maskmpl = false(npxt,nl);      % Posterior lateral mask
      end
%
      idapm = sum(iap(~idl,:));
      idam = idapm<0;                  % Column index to anterior medial masks
      if any(idam)
        maskmam = maskmm(:,:,idam);    % Anterior medial mask
      else
        maskmam = false(npxt,nm);      % Anterior medial mask
      end
      idpm = idapm>0;                  % Column index to posterior medial masks
      if any(idpm)
        maskmpm = maskmm(:,:,idpm);    % Posterior medial mask
      else
        maskmpm = false(npxt,nm);      % Posterior medial mask
      end
%
% Get Compartment and AP Specific Meniscus Coordinates
%
      mcoordl = mcoords(idl,:);        % Meniscus coordinates for lateral slices
      mcoordm = mcoords(~idl,:);       % Meniscus coordinates for medial slices
%
      mcoordal = mcoordl(:,idal);      % Anterior lateral meniscus coordinates
      mcoordpl = mcoordl(:,idpl);      % Posterior lateral meniscus coordinates
%
      mcoordam = mcoordm(:,idam);      % Anterior medial meniscus coordinates
      mcoordpm = mcoordm(:,idpm);      % Posterior medial meniscus coordinates
%
% Save Masks, ROIS and Slice Information into MAT File
%
      savnam = [mnam(1:end-4) '_mrois.mat'];
      savnam = fullfile(sdir,savnam);
      save(savnam,'bbox','cmprt','legs','maskmal','maskmam', ...
           'maskmpl','maskmpm','masks','mcoordal','mcoordam', ...
           'mcoordpl','mcoordpm','rsl','rsll','rslm');
%
      close all;
%
   end                  % End of m loop - MAT file loop
%
% end                     % End of subject (ks) loop
%
return