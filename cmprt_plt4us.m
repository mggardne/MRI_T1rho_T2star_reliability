function rimgr = cmprt_plt4us(v,mask1,rsls1,mask2,rsls2,idt,slk, ...
                              tcp1,nps1,tcp2,nps2,mxtc)
%CMPRT_PLT4US  Creates on image matrix with T1rho/T2* values for the
%          underlying images for a single slice within two regions of
%          interest (ROIs).
%
%          RIMGR = CMPRT_PLT4U(V,MASK1,RSLS1,MASK2,RSLS2,IDT,SLK,TCP1,
%          NPS1,TCP2,NPS2) Given a four-dimensional matrix of T1/T2
%          intensities from a MRI image volume, V, where the first two
%          dimensions are an image, the third dimension are the slices,
%          and the fourth dimension are the spin lock/echo times, three
%          dimensional logical masks with the first dimension being the
%          image, the second dimension being the superficial layer in
%          the first column and deep layer in the second column and the
%          third dimension being slices in a cell array of masks for the
%          first regions of interest (ROIs) with the first index to the
%          lateral and medial compartments and the second index to the
%          femur and tibia, MASK1, a cell array with the slices within
%          each compartment, RSLS1, a cell array of masks for the
%          second ROIs with an index to the lateral and medial
%          compartments, MASK2, a cell array with the slices within
%          each compartment, RSLS2, index to the spin lock/echo time to
%          use for plotting, IDT, cell array of T1rho/T2* values for the
%          first ROIs, SLK, slice number to be plotted, TCP1, the number
%          of fitted pixels in each slice within the compartments in a
%          cell array, NPS1, cell array of T1rho/T2* values for the
%          second ROIs, TCP2, and the number of fitted pixels in each
%          slice within the compartments in a cell array, NPS2, plots
%          the T1rho/T2* values with the underlying images by slice
%          within the regions of interest defined by the MASK1 and
%          MASK2.
%
%          CMPRT_PLT4U(V,MASK1,RSLS1,MASK2,RSLS2,IDT,SLK,TCP1,NPS1,TCP2,
%          NPS2,MXTC) Uses MXTC as the maximum plotting value for the
%          color.
%
%          NOTES:  1.  Plots the pixel output of cmprt_ana4.m and
%                  cmprt_ana4u.m.  See cmprt_ana4.m, cmprt_ana4u.m, and
%                  mri_fitr4u.m.
%
%          05-Dec-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<11
  error([' *** ERROR in cmprt_plt4us:  Eleven input variables are', ...
         ' required!']);
end
%
if nargin<11||isempty(mxtc)
  mxtc = 70;
end
%
% Initialize Arrays
%
idxs1 = [2 2];          % Maximum indices for compartment and bone
%
% Find Compartment
%
if any(rsls1{1}==slk)
  kr = 1;               % Lateral compartment
else
  kr = 2;               % Medial compartment
end
%
% Get Compartment Variables
%
mskr1 = mask1{kr};      % Region 1 mask for this compartment
mskr2 = mask2{kr};      % Region 2 mask for this compartment
rsl1 = rsls1{kr};       % Region 1 slices for this compartment
rsl2 = rsls2{kr};       % Region 2 slices for this compartment
ks1 = find(rsl1==slk);  % Slice in region 1
ks2 = find(rsl2==slk);  % Slice in region 2
%
% Get Slice Image
%
rimg = squeeze(v(:,:,slk,idt));   % T1/T2 data for slice and plot spin lock/echo time
rimgr = rimg;           % Image for results
%
% Scale T1rho/T2* Image to -mxtc to Zero
%      
rimgr = rimgr-min(rimgr(:));
imgmx = max(rimgr);
rimgr = mxtc*rimgr./imgmx;             % Scale MRI image
rimgr = rimgr-(mxtc+0.01);
%
% Get Region 1 Results
% Loop through Bone
%
for kb = 1:2
%
   mskb1 = mskr1{kb};   % Mask for this bone
%
   msk1 = squeeze(mskb1(:,:,ks1));          % Mask for this slice
   msk1 = squeeze(msk1(:,1)|msk1(:,2));     % Combine cartilage layers
%
   idx1 = sub2ind(idxs1,kb,kr);        % Index to T1rho/T2* results
   npsk1 = nps1{idx1};
   npsks1 = sum(npsk1(1:ks1));
   npsks1 = (npsks1-npsk1(ks1)+1:npsks1)';
%
   tcpk1 = tcp1{idx1};
%
   rimgr(msk1) = tcpk1(npsks1);        % T1rho/T2* values
%
end             % End of kb loop - bones loop
%
% Get Region 2 Results
%
msk2 = squeeze(mskr2(:,:,ks2));        % Mask for this slice
msk2 = squeeze(msk2(:,1)|msk2(:,2));   % Combine cartilage layers
%
npsk2 = nps2{kr};
npsks2 = sum(npsk2(1:ks2));
npsks2 = (npsks2-npsk2(ks2)+1:npsks2)';
%
tcpk2 = tcp2{kr};
%
rimgr(msk2) = tcpk2(npsks2);           % T1rho/T2* values
%
return