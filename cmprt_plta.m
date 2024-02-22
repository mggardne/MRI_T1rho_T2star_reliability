function cmprt_plta(v,mask,rsls,nrsls,idt,tcp,nps,mxtc,cmap,txt1,psnam)
%CMPRT_PLTA Plots T1rho/T2* values with the underlying images by slice 
%          within regions of interest (ROIs) defined by a logical mask.
%
%          CMPRT_PLTA(V,MASK,RSLS,NRSLS,IDT,TCP,NPS) Given a four-
%          dimensional matrix of T1/T2 intensities from a MRI image
%          volume, V, where the first two dimensions are an image, the
%          third dimension are the slices and the fourth dimension are
%          the spin lock/echo times, three dimensional logical mask
%          with the first dimension being the image, the second
%          dimension being both layers combined in the first column, the
%          superficial layer in the second column and deep layer in the
%          third column, and the third dimension being slices, MASK, a
%          cell array with the slices within each compartment, RSLS, the
%          number of slices in each compartment, NRSLS, index to the
%          spin lock/echo time to use for plotting, IDT, cell array of
%          T1rho/T2* values, TCP, and the number of fitted pixels in
%          each slice within the compartments in a cell array, NPS,
%          plots the T1rho/T2* values with the underlying images by
%          slice within the regions of interest defined by the MASK.
%
%          CMPRT_PLTA(V,MASK,RSLS,NRSLS,IDT,TCP,NPS,MXTC,CMAP,TXT1)
%          Given the maximum plotting value for the color scale, MXTC, a
%          three column color map, CMAP, and a text string for the first
%          line of the plot title, TXT1, plots the T1rho/T2* values with
%          a color maximum of MXTC using the color map, CMAP, and using
%          TXT1 for the first line of the plot title.  The default
%          maximum value is 70.  The default color map is gray for the
%          image and jet for the T1rho/T2* values.  The default first
%          line title text is "Results Plot".
%
%          CMPRT_PLTA(V,MASK,RSLS,NRSLS,IDT,TCP,NPS,MXTC,CMAP,TXT1,
%          PSNAM) Given the name for a PS file, PSNAM, prints the plots
%          to the PS file.  By default, the plots are not printed.
%
%          NOTES:  1.  Plots the pixel output of cmprt_ana_all.m.  See
%                  cmprt_ana_all.m and mri_fitra.m.
%
%          09-Nov-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<7
  error([' *** ERROR in cmprt_plta:  Seven input variables are', ...
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
if nargin<10||isempty(txt1)
  txt1 = 'Results Plot';
end
%
if nargin<11||isempty(psnam)
  isave = false;
else
  isave = true;
end
%
% Compartment Labels
%
tcmprts = {'Lateral'; 'Medial'};
%
% Get Combined Layers Mask
%
mskr = squeeze(mask(:,1,:));           % Mask for all slices
%
rsl = rsls{3};          % Slices for all compartments
nrsl = nrsls(3);        % Number of slices
rslm = rsls{2};         % Medial compartment slices
%
idxs = [3 3 3];         % Maximum indices for compartment, bone and layer
%
% Loop through Slices
%
for ks = 1:nrsl
%
   msk = mskr(:,ks);    % Mask for this slice
%
   slk = rsl(ks);       % Slice
%
% Get Slice Image
%
   rimgr = squeeze(v(:,:,slk,idt));    % T1/T2 data for slice and plot spin lock/echo time
%
% Scale T1rho/T2* Image to -mxtc to Zero
%      
   rimgr = rimgr-min(rimgr(:));
   imgmx = max(rimgr);
   rimgr = mxtc*rimgr./imgmx;
   rimgr = rimgr-(mxtc+0.01);
%
   npsks = sum(nps(1:ks));
   npsks = (npsks-nps(ks)+1:npsks)';
%
   rimgr(msk) = tcp(npsks);            % T1rho/T2* values
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
%
   idc = ismember(slk,rslm)+1;
   tcmprt = [tcmprts{idc} ' Compartment'];
%
   if iscell(txt1)
     ttxt = {txt1{:} [tcmprt ', Slice ' int2str(slk)]}';
   else
     ttxt = {txt1; [tcmprt ', Slice ' int2str(slk)]};
   end
%
   title(ttxt,'FontSize',16,'FontWeight','bold');
%
   hb = colorbar;
   set(hb,'Limits',[0 mxtc]);
%
   if isave             % Print plots
     if ks==1
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
return