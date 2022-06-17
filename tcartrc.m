%#######################################################################
%
%       * Tibia CARilage Thickness Reliability Combined Program *
%
%          M-File which calculates the tibia cartilage thicknesses for
%     the 10 reliability subjects' second visit unloaded left knee using
%     both the T1 FFE and T1rho MRI images.  The program reads the tibia
%     CVS files from the subject visit 2 subdirectories and calculates
%     the origin and rotation matrix to transform the MRI data to a
%     tibia coordinate system, defines a 1 mm by 1 mm grid, scales and
%     projects the grid onto the bone surface mesh, calculates the
%     intersection of the bone normals with the cartilage surface mesh
%     and calculates cartilage thicknesses.  The same grid is projected
%     onto both MRI acquisitions (T1 FFE and T1rho) so similar points
%     can be compared across the two acquisitions.
%
%     NOTES:
%
%     1.  The M-files car_thk8.m, comp_msh.m, coord_tf.m,
%     get_subj.m, gridproj.m, in_tri2d.m, inert_tri.m, isect.m,
%     line_fit.m, mk_tri4a.m, mk_tri4i.m, mk_tri4p2.m,
%     mk_tri4s.m, nod2tri.m, nod_norm.m, pcsr.m, plane_fit.m,
%     pts2lin.m, rd_roi3.m, sl_info.m, tri_fix2.m, tri_norm.m,
%     tsect4.m and xprod.m must be in the current path or
%     directory.
%
%     2.  Subject directories must be two character numbers (i.e. 01,
%     02, 05, ..., 10, 12, 13).  No other directories should start with
%     a "0" or "1".
%
%     3.  CSV data files must be in the subject subdirectories:
%     "Visit 2\RHO\TIBIA" or "Visit 2\FFE".
%
%     4.  CSV data file names must be in a specific format.
%     The format is:  subject number _L_ {AX or SAG (for bone) or
%     SAGAR (cartilage)} _TIB_ {FFE or T1R} _V2_*.csv
%
%     5.  Results are put in the directory "TibCartThk".  This program
%     generates MAT files with bone coordinates and bone meshes with
%     file names:  [subject IDs _ {FFE or T1R} '_bone.mat'].  Produces
%     a MAT file with the grid and scaling for the tibias:  tgridrc.mat.
%     Also, produces MAT files with cartilage coordinates and cartilage
%     meshes with file names:  [subject IDs _ {FFE or T1R} '_cart.mat'].
%     Final thicknesses are saved in the MAT file:  tcartrcs.mat.
%
%     11-May-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Plot Parameter
%
% iplt = false;   % No plots
iplt = true;    % Plots
%
% iprt = false;   % No printing of plots
iprt = true;    % Print plots
%
% Get Two Character Subject Directories
%
rdirs = dir;
rdirs = {rdirs([rdirs.isdir]').name}';
idd = startsWith(rdirs,'0')|startsWith(rdirs,'1');
rdirs = rdirs(idd);
idv = cellfun('length',rdirs)==2;      % Two character subject directories
rdirs = char(rdirs(idv));
ntib = size(rdirs,1);   % Number of subjects/tibias
%
% Initialize Coordinate and Bone Arrays
%
leg = 'L';              % Only left legs
lnam = 'Left Leg';      % Leg name
legn = 0;               % Leg number (0 = left)
%
xyz0 = zeros(ntib,3,2); % Coordinate system origins
rotmat = zeros(3,3,ntib,2);       % Rotation matrices (MRI to PCS)
tribs = cell(ntib,2);   % Sagittal bone triangular mesh connectivity
xyzbs = cell(ntib,2);   % Sagittal bone point coordinates
xyzmns = zeros(ntib,3,2);         % Minimum individual tibia coordinates
xyzmxs = zeros(ntib,3,2);         % Maximum individual tibia coordinates
%
% Loop through Subjects/Tibias
%
for kt = 1:ntib
%
% Get Subject Number
%
   rdir = rdirs(kt,:);
   subjn = eval(rdir);  % Subject number
%
   bfnam = [fnam '_bone.mat']; % Bone MAT file name
%
% Get or Generate Bone Data
%
   if exist(bfnam,'file')
     load(bfnam);
     cgs(k,:) = cg;
     rotmat(:,:,k) = r;
     triss{k} = tris;
     xyzss{k} = xyzs;
   else
%
% Get Bone File Names
%
     d = dir(fullfile(csvdir,[fnam '_AX_PATBONE_*.csv']));      % Axial file name
     fnam_ax = d.name;  % Axial file name
%
     d = dir(fullfile(csvdir,[fnam '_SAG_PATBONE_*.csv']));     % Sagittal file name
     fnam_sag = d.name; % Sagittal file name
%
% Read Bone Files and Get Data Cell Arrays
%
     roiax = rd_roi3(fullfile(csvdir,fnam_ax));
     roisag = rd_roi3(fullfile(csvdir,fnam_sag));
%
     datx = roiax.data';       % Axial data
     dats = roisag.data';      % Sagittal data
%
% Get Patella Coordinate System
%
     [cgs(k,:),rotmat(:,:,k),xyzrs(k,:),vecrs(k,:),vols(k), ...
      triss{k},xyzss{k}] = pcsr(kleg,dats,datx,iplt);
%
     if iplt
       hf = get(0,'Children');
       fn = cell2mat(get(hf,'Number'));
       [~,ids] = sort(fn);
       hf = hf(ids);    % Figure handles in order
%
       if iprt
         psnam = ['pcsrc3_' kid '_r.ps'];
       end
       hca = get(hf,'child');  % Handles to plot axes
       for l = 1:length(hca)
          ht = get(hca{l},'Title');    % Get title handles
          htt = get(ht,'String');
          set(ht,'String',{['Knee ' kid]; ['CS R ' htt]}, ...
              'Interpreter','none');   % Update titles
          if iprt
            if l==1
              print(hf(l),'-dpsc2',psnam);
            else
              print(hf(l),'-dpsc2','-append',psnam);
            end
          end
       end
%
       close all;
%
     end