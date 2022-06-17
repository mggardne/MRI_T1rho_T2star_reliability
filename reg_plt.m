function reg_plt(v,mask,rsls,nrsls,idt,tcp,nps,mxtc,cmap,txt1,psnam)
%REG_PLT   Plots T1rho/T2* values with the underlaying images by slice 
%          within regions of interest.
%
%          REG_PLT(V,MASK,RSLS,NRSLS,IDT,TCP,NPS) Given a four-
%          dimensional matrix of T1/T2 intensities from a MRI image
%          volume, V, where the first two dimensions are an image, the
%          third dimension are the slices, and the fourth dimension are
%          the spin lock/echo times, three dimensional logical masks
%          with the first dimension being the image, the second
%          dimension being the superficial layer in the first column and
%          deep layer in the second column and the third dimension being
%          slices in a cell array of masks with the first index to
%          regions 1 and 2 and the second index to the femur and tibia,
%          MASK, a cell array with the slices within each region, RSLS,
%          the number of slices in each region, NRSLS, index to the
%          spin lock/echo time to use for plotting, IDT, cell array of
%          T1rho/T2* values, TCP, and the number of fitted pixels in
%          each slice within the regions in a cell array, NPS, plots
%          the T1rho/T2* values with the underlaying images by slice 
%          within the regions of interest defined by the MASK.
%
%          REG_PLT(V,MASK,RSLS,NRSLS,IDT,TCP,NPS,MXTC,CMAP,TXT1) Given
%          the maximum plotting value for the color scale, MXTC, a three
%          color map, CMAP, and a text string for the first line of the
%          plot title, TXT1, plots the T1rho/T2* values with a color
%          maximum of MXTC using the color map, CMAP, and using TXT1 for
%          the first line of the plot title.  The default maximum value
%          is 70.  The default color map is gray for the image and jet
%          for the T1rho/T2* values.  The default first line title text
%          is "Results Plot".
%
%          REG_PLT(V,MASK,RSLS,NRSLS,IDT,TCP,NPS,MXTC,CMAP,TXT1,PSNAM)
%          Given the name for a PS file, PSNAM, prints the plots to the
%          PS file.  By default, the plots are not printed.
%
%          NOTES:  1.  Plots the pixel output of reg_ana.m.  See
%                  reg_ana.m and mri_fitr.m.
%
%          26-Jan-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<7
  error([' *** ERROR in reg_plt:  Seven input variables are', ...
         ' required!']);
end
%
if nargin<8||isempty(mxtc)
  mxtc = 70;
end
%
if nargin<9||isempty(cmap)
%
% Default Color Map
%
  gmap = gray(128);     % Gray color map for not cartilage
  jmap = jet(128);      % Jet color map for cartilage measures
  cmap = [gmap; jmap];
end
%
if nargin<1||isempty(txt1)
  txt1 = 'Results Plot';
end
%
if nargin<11||isempty(psnam)
  isave = false;
else
  isave = true;
end
%
% Initialize Arrays
%
idxs = [2 2 2];         % Maximum indices for region, bone and layer
%
% Loop through Regions
%
for kr = 1:2
%
   mskr = mask{kr};     % Mask for this region
   rsl = rsls{kr};      % Slices for this region
   nrsl = nrsls(kr);    % Number of slices in this region
%
% Loop through Slices
%
   for ks = 1:nrsl
%
      slk = rsl(ks);    % Slice
%
% Get Slice Image
%
      rimg = squeeze(v(:,:,slk,idt));   % T1/T2 data for slice and plot spin lock/echo time
      rimgp = rimg;     % Image for plotting regions
      rimgr = rimg;     % Image for plotting results
%
% Scale T2 Image to -mxtc to Zero
%      
      rimgr = rimgr-min(rimgr(:));
      imgmx = max(rimgr);
      rimgr = mxtc*rimgr./imgmx;
      rimgr = rimgr-(mxtc+0.01);
%
% Loop through Bone
%
      for kb = 1:2
%
         mskb = mskr{kb};              % Mask for this bone
%
% Loop through Layer
%
         for kl = 1:2
%
            msk = squeeze(mskb(:,kl,ks));   % Mask for this slice and layer
%
            idx = sub2ind(idxs,kl,kb,kr);   % Index to T1rho/T2* results
            npsk = nps{idx};
            npsks = sum(npsk(1:ks));
            npsks = (npsks-npsk(ks)+1:npsks)';
%
            tcpk = tcp{idx};
%
            rimgr(msk) = tcpk(npsks);  % T1rho/T2* values
%
         end            % End of kl loop - layers loop
      end               % End of kb loop - bones loop
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
      title({txt1; ['Slice ' int2str(slk)]},'FontSize',16, ...
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
   close all;
%
end                     % End of kr loop - regions loop
%
return