function cmprt_pltu(v,mask1,rsls1,mask2,rsls2,idt,tcp1,nps1,tcp2, ...
                    nps2,mxtc,cmap,txt1,psnam)
%CMPRT_PLTU  Plots T1rho/T2* values with the underlying images by slice 
%          within regions of interest (ROIs).
%
%          CMPRT_PLTU(V,MASK1,RSLS1,MASK2,RSLS2,IDT,TCP1,NPS1,TCP2,NPS2)
%          Given a four-dimensional matrix of T1/T2 intensities from
%          a MRI image volume, V, where the first two dimensions are an
%          image, the third dimension are the slices, and the fourth
%          dimension are the spin lock/echo times, three dimensional
%          logical masks with the first dimension being the image, the
%          second dimension being the superficial layer in the first
%          column and deep layer in the second column and the third
%          dimension being slices in a cell array of masks for the
%          first regions of interest (ROIs) with the first index to the
%          lateral and medial compartments and the second index to the
%          femur and tibia, MASK1, a cell array with the slices within
%          each compartment, RSLS1, a cell array of masks for the
%          second ROIs with an index to the lateral and medial
%          compartments, MASK2, a cell array with the slices within
%          each compartment, RSLS2, index to the spin lock/echo time to
%          use for plotting, IDT, cell array of T1rho/T2* values for the
%          first ROIs, TCP1, the number of fitted pixels in each slice
%          within the compartments in a cell array, NPS1, cell array of
%          T1rho/T2* values for the second ROIs, TCP2, and the number
%          of fitted pixels in each slice within the compartments in a
%          cell array, NPS2, plots the T1rho/T2* values with the
%          underlying images by slice within the regions of interest
%          defined by the MASK1 and MASK2.
%
%          CMPRT_PLTU(V,MASK1,RSLS1,MASK2,RSLS2,IDT,TCP1,NPS1,TCP2,NPS2,
%          MXTC,CMAP,TXT1) Given the maximum plotting value for the
%          color scale, MXTC, a three color map, CMAP, and a text
%          string for the first line of the plot title, TXT1, plots the
%          T1rho/T2* values with a color maximum of MXTC using the
%          color map, CMAP, and using TXT1 for the first line of the
%          plot title.  The default maximum value is 70. The default
%          color map is gray for the image and jet for the T1rho/T2*
%          values.  The default first line title text is "Results Plot".
%
%          CMPRT_PLTU(V,MASK1,RSLS1,MASK2,RSLS2,IDT,TCP1,NPS1,TCP2,NPS2,
%          MXTC,CMAP,TXT1,PSNAM)  Given the name for a Postscript (PS)
%          file, PSNAM, prints the plots to the PS file.  By default,
%          the plots are not printed.
%
%          NOTES:  1.  Plots the pixel output of cmprt_ana.m and
%                  cmprt_anau.m.  See cmprt_ana.m, cmprt_anau.m, and
%                  mri_fitr3u.m.
%
%          21-Oct-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<10
  error([' *** ERROR in cmprt_pltu:  Ten input variables are', ...
         ' required!']);
end
%
if nargin<11||isempty(mxtc)
  mxtc = 70;
end
%
if nargin<12||isempty(cmap)
%
% Default Color Map
%
  gmap = gray(128);     % Gray color map for not cartilage
  jmap = jet(128);      % Jet color map for cartilage measures
  cmap = [gmap; jmap];
end
%
if nargin<13||isempty(txt1)
  txt1 = 'Results Plot';
end
%
if nargin<14||isempty(psnam)
  isave = false;
else
  isave = true;
end
%
% Compartment Labels
%
tcmprts = {'Lateral'; 'Medial'};
%
% Initialize Arrays
%
idxs1 = [2 2 2];        % Maximum indices for compartment, bone, and layer
idxs2 = [2 2];          % Maximum indices for compartment and layer
%
% Loop through Compartments
%
for kr = 1:2
%
   mskr1 = mask1{kr};   % Region 1 mask for this compartment
   mskr2 = mask2{kr};   % Region 2 mask for this compartment
   rsl1 = rsls1{kr};    % Region 1 slices for this compartment
   rsl2 = rsls2{kr};    % Region 2 slices for this compartment
   rslp = unique([rsl1; rsl2]);
   nrslp = size(rslp,1);     % Number of slices in this compartment
   tcmprt = [tcmprts{kr} ' Compartment'];
%
% Loop through Slices
%
   ks1 = 0;             % Index to Region 1 slices
   ks2 = 0;             % Index to Region 2 slices
%
   for ks = 1:nrslp
%
      slk = rslp(ks);   % Slice
%
% Get Slice Image
%
      rimg = squeeze(v(:,:,slk,idt));   % T1/T2 data for slice and plot spin lock/echo time
      rimgr = rimg;     % Image for plotting results
%
% Scale T1rho/T2* Image to -mxtc to Zero
%      
      rimgr = rimgr-min(rimgr(:));
      imgmx = max(rimgr);
      rimgr = mxtc*rimgr./imgmx;
      rimgr = rimgr-(mxtc+0.01);
%
% Get Region 1 Results
%
      if any(rsl1==slk)
%
        ks1 = ks1+1;    % Region 1 slices
%
% Loop through Bone
%
        for kb = 1:2
%
           mskb1 = mskr1{kb};          % Mask for this bone
%
% Loop through Layer
%
           for kl = 1:2 % 1 - superficial and 2 - deep
%
              msk1 = squeeze(mskb1(:,kl,ks1));   % Mask for this slice and layer
%
              idx1 = sub2ind(idxs1,kl,kb,kr);    % Index to T1rho/T2* results
              npsk1 = nps1{idx1};
              npsks1 = sum(npsk1(1:ks1));
              npsks1 = (npsks1-npsk1(ks1)+1:npsks1)';
%
              tcpk1 = tcp1{idx1};
%
              rimgr(msk1) = tcpk1(npsks1);       % T1rho/T2* values
%
           end          % End of kl loop - layers loop
        end             % End of kb loop - bones loop
      end               % End of if Region 1
%
% Get Region 2 Results
%
      if any(rsl2==slk)
%
        ks2 = ks2+1;    % Region 2 slices
%
% Loop through Layer
%
        for kl = 1:2    % 1 - superficial and 2 - deep
%
           msk2 = squeeze(mskr2(:,kl,ks2));      % Mask for this slice and layer
%
           idx2 = sub2ind(idxs2,kl,kr);          % Index to T1rho/T2* results
           npsk2 = nps2{idx2};
           npsks2 = sum(npsk2(1:ks2));
           npsks2 = (npsks2-npsk2(ks2)+1:npsks2)';
%
           tcpk2 = tcp2{idx2};
%
           rimgr(msk2) = tcpk2(npsks2);          % T1rho/T2* values
%
        end             % End of kl loop - layers loop
      end               % End of if Region 2
%
% Plot Slice
%
      figure;
      orient landscape;
%
      imagesc(rimgr,[-mxtc mxtc]);
      colormap(cmap);
      axis image;
      axis off;
      title({txt1; ['Slice ' int2str(slk)]; tcmprt},'FontSize',16, ...
            'FontWeight','bold');
%
      hb = colorbar;
      set(hb,'Limits',[0 mxtc]);
%
      if isave          % Print plots
        if kr==1&&ks==1
          print('-dpsc2','-r600','-fillpage',psnam);
        else
          print('-dpsc2','-r600','-fillpage','-append',psnam);
        end
      end
%
   end                  % End of ks loop - slices loop
%
%    close all;
%
end                     % End of kr loop - compartments loop
%
return