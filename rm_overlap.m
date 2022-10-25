%#######################################################################
%
%             * ReMove OVERLAPs in Segmentations Program *
%
%          M-File which reads the MS-Excel spreadsheet file,
%     contact_chk.xlsx, to find the segmentation CSV files and slices
%     where the femur and tibia segmentations overlap, reads the femur
%     and tibia segmentation CSV files, finds and adds midpoints
%     between the femur and tibia segmentations and writes the revised
%     data to new segmentation CSV files.  The new segmentation CSV
%     files have the same name with "_RO" appended to the name.
%
%          The slices with corrected overlaps are plotted to Postscript
%     files rm_overlap*.ps in the Results\RemoveOverlap folder.
%
%     NOTES:  1.  M-files decomp.m, dis2lin.m, get_csv.m, get_olap.m,
%             lisect3.m, lisect4a.m, lisect5.m, lsect3.m, lsect4.m,
%             lsect5.m, mk_contct.m, mxd2lins.m, parse_csv.m, pts2lin.m,
%             rd_roi6.m, rd_rois.m, trnsf2mm.m and trnsf2pixel.m must be
%             in the current directory or path.
%
%             2.  Only files with overlaps greater than a minimum
%             overlap are corrected.  See variable "omin" on line 37.
%
%             3.  For a couple of slices, the 2D and 3D intersections
%             do not agree.  Future versions should compare the 2D and
%             3D points and verify that the 2D and 3D corrections are
%             the same.  If not, the user should pick either the 2D or
%             3D correction based on the plots.
%
%     25-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Analysis Parameters
%
omin = 0.2;             % Minimum overlap for correction in pixels
%
tol2 = 1e-4;            % Tolerance for 2D intersections
tol3 = 2.5e-1;          % Tolerance for 3D intersections
%
% Print?
%
% iprt = false;
iprt = true;
%
% Output Directory and Postscript Filename
%
resdir = fullfile('Results','RemoveOverlap');    % Results directory
psnam = fullfile(resdir,'rm_overlap_');     % Start of PS file name
pstyp = '.ps';          % PS file type
%
% Read Overlap Data from contact_chk.xlsx
%
% csvdat
%     Column 1:  Subject
%     Column 2:  Visit
%     Column 3:  Result
%     Column 4:  Leg
%     Column 5:  Load
%
datdir = fullfile('Results','ContactChecks');
xlsnam = 'contact_chk.xlsx';
% sht = 'Combined Data';
sht = 1;
[odat,csvdat,~,sls,ids] = get_olap(xlsnam,sht,datdir);
%
% Get CSV Files with Overlaps Greater than the Minimum Overlap
%
ovlap = odat(:,8);      % Femur-tibia overlaps in pixels
idcorr = ovlap>omin;    % Logical index to overlaps to be corrected (> minimum)
%
sls = sls(idcorr,:);    % Slices with overlaps greater than the minimum
ids = ids(idcorr,:);    % Index to slices within files
sids = sum(ids);
idcorr = sids>0;        % Index to files with overlaps > the minimum 
ids = ids(:,idcorr);    % Update slices index
%
csvdat = csvdat(idcorr,:);             % Update files with overlaps
ncsv = size(csvdat,1);  % Number of files with overlaps to be corrected
%
% Get Subject Directories, Subject Numbers and Visit Subdirectories
%
sdirs = dir('MRIR*');
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
%
stxts = char(sdirs);
stxts = stxts(:,6:7);   % Subject numbers as text
subj_nums = str2double(cellstr(stxts)); % Subject numbers
%
vdirs = {'Visit1'; 'Visit2'};          % Visit directories
vnams = ['Visit 1'; 'Visit 2'];
%
% Get Analysis and Bone Subdirectories
%
adirs = ['RHO'; 'T2S'];                % Analysis directories
anams = {'T1\rho'; 'T2*'};
%
irhos = [4; 5];         % Number of spin lock (4) or echo (5) times
%
bdirs = ['Femur'; 'Tibia'];            % Bone directories
%
% Get Legs and Loads Variables
%
legs = ['L'; 'R'];      % Left/Right
lnams = {'Left Leg'; 'Right Leg'};
%
lds = ['UL'; 'LD'];     % Unloaded/loaded
ldnams = {'Unloaded'; 'Loaded'};
%
% Loop through CSV Files
%
% for ko = 1:ncsv
% for ko = 27:ncsv        % Start with subject 5
% for ko = 100:ncsv       % Start with subject 10
% for ko = 113:ncsv       % Start with subject 10's 14th CSV file
% for ko = 119:119        % Do subject 12's 5th CSV file
% for ko = 67:67          % Do subject 7's 13th CSV file
for ko = 111:114        % Do subject 10's 2nd Visit T2* CSV files
%
% Get CSV Identifiers and Parse the Identifiers
%
   csvk = csvdat(ko,:); % CSV file identifiers
%
   idsn = subj_nums==csvk(1);          % Get index to subject number
   sdir = sdirs{idsn};                 % Subject directory
   stxt = stxts(idsn,:);               % Subject number as text
%
   psnams = [psnam stxt];              % Add subject to PS file name
%
   leg = legs(csvk(4)); % Leg
   ltxt = lnams{csvk(4)};
%
   ld = lds(csvk(5),:); % Load
   ldtxt = ldnams{csvk(5)};
%
% Get Slices
%
   slsk = sls(ids(:,ko));              % Slices for this subject
   ns = size(slsk,1);   % Number of slices
%
% Get Directories
%
   cdir = fullfile(sdir,vdirs{csvk(2)},adirs(csvk(3),:));
%
   vtxt = vnams(csvk(2),:);
   atxt = anams{csvk(3)};
%
   psnamv = [psnams '_V' int2str(csvk(2))]; % Add visit to PS file name
   psnama = [psnamv '_' adirs(csvk(3),:)];  % Add result type to PS file name
   psnamf = [psnama '_' leg '_' ld pstyp];  % Add leg and load to PS file name
%
   cdirf = fullfile(cdir,bdirs(1,:));  % Femur segmentation CSV directory
   cdirt = fullfile(cdir,bdirs(2,:));  % Tibia segmentation CSV directory
%
% Get CSV Files
%
   csvf = get_csv(cdirf,leg,ld);       % Femur CSV file
   csvt = get_csv(cdirt,leg,ld);       % Tibia CSV file
%
% Read CSV File and Get Lines with Slices
%
   irho = irhos(csvk(3));
%
   [flines,fhdrs,fids] = rd_csv(cdirf,csvf,slsk,irho);     % Femur
   [tlines,thdrs,tids] = rd_csv(cdirt,csvt,slsk,irho);     % Tibia
%
% Loop through Slices
%
   for ks = 1:ns
%
% Get Segmentation Points from Slice Lines
%
      [fpts3d,fpts2d,~,fidn,fidx] = parse_csv(fhdrs, ...
                                              flines{fids(ks)});     % Femur
      [tpts3d,tpts2d,~,tidn,tidx] = parse_csv(thdrs, ...
                                              tlines{tids(ks)});     % Tibia
%
% Find Intersections Between the Tibia and Femur Segmentations
%
      [ipts2d,~,idcs2d] = lsect5(tpts2d,fpts2d,tol2);
      [ipts3d,~,idcs3d] = lisect5(tpts3d,fpts3d,tol3);
%
% Make Sure Intersections are in Order
%
      [ipts2d,idsrt] = sortrows(ipts2d,1);
      idcs2d = idcs2d(idsrt,:);
%
      [ipts3d,idsrt] = sortrows(ipts3d,2);
      idcs3d = idcs3d(idsrt,:);
%
% Remove the Overlaps in the 2D and 3D Point Data
%
      line1 = ['Subject ' stxt ' ' vtxt ' ' atxt];
      line2 = [ltxt ', ' ldtxt ', Slice ' int2str(slsk(ks))];
      ttxt = {[line1 ' 2D']; line2};   % Plot title text
%
      [tpts2dro,fpts2dro] = mk_contct(tpts2d,fpts2d,ipts2d,idcs2d, ...
                                      true,ttxt);
      if iprt
        if ks==1
          print('-dpsc2','-r600','-fillpage',psnamf);
        else
          print('-dpsc2','-r600','-fillpage','-append',psnamf);
        end
      end
%
      ttxt = {[line1 ' 3D']; line2};   % Plot title text
%
      [tpts3dro,fpts3dro] = mk_contct(tpts3d,fpts3d,ipts3d,idcs3d, ...
                                      true,ttxt);
      if iprt
        print('-dpsc2','-r600','-fillpage','-append',psnamf);
      end
%
% Update Segmentation File Lines with the New 3D and 2D Points
%
      npts = size(fpts3dro,1);
      npts2d = size(fpts2dro,1);
      if npts~=npts2d
        warning([' *** WARNING in RM_OVERLAP:  Numbers of 3D and', ...
                 ' 2D points do not match!']);
        if ko==119      % Based on plots
          fprintf(1,['  Using transformed 2D points for 3D points', ...
                     ' for CSV file: %s\n\n'],csvf);
          fpts3dro = trnsf2mm(fpts2d,fpts3d,fpts2dro);
          fprintf(1,['  Using transformed 2D points for 3D points', ...
                     ' for CSV file: %s\n\n'],csvt);
          tpts3dro = trnsf2mm(tpts2d,tpts3d,tpts2dro);
        else
          fprintf(1,['  Using transformed 3D points for 2D points', ...
                     ' for CSV file: %s\n\n'],csvf);
          fpts2dro = trnsf2pixel(fpts3d,fpts2d,fpts3dro);
        end
      end
%
      flines{fids(ks)} = line_upd(flines{fids(ks)},fidn,fidx, ...
                                  fpts3dro,fpts2dro);
%
      npts = size(tpts3dro,1);
      npts2d = size(tpts2dro,1);
      if npts~=npts2d
        warning([' *** WARNING in RM_OVERLAP:  Numbers of 3D and', ...
                 ' 2D points do not match!']);
        fprintf(1,['  Using transformed 3D points for 2D points', ...
                   ' for CSV file: %s\n\n'],csvt);
        tpts2dro = trnsf2pixel(tpts3d,tpts2d,tpts3dro);
      end
%
      tlines{tids(ks)} = line_upd(tlines{tids(ks)},tidn,tidx, ...
                                  tpts3dro,tpts2dro);
%
%       pause;
%       pause(0.1);
      close all;
%
   end
%
% Write Updated Segmentation lines to a New CSV File 
%
   idot = strfind(csvf,'.');
   idot = idot(end);
   csvf = [csvf(1:idot-1) '_RO' csvf(idot:end)];
   write_csv(cdirf,csvf,flines);
%
   idot = strfind(csvt,'.');
   idot = idot(end);
   csvt = [csvt(1:idot-1) '_RO' csvt(idot:end)];
   write_csv(cdirt,csvt,tlines);
%
end
%
return