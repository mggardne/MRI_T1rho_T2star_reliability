%#######################################################################
%
% * FIX Subject 7 Remove Overlap (RO) 2 Segmentation CSV File Program *
%
%          M-File which reads and fixes a segmentation CSV file where
%     the 3D and 2D coordinates in the femur CSV segmentation file are
%     not similar.  The program only checks the single CSV files with
%     "_RO2" in the filename (created by rm_overlapMay24.m):
%     "007_L_SAGAR_FEM_T2S_LD_V2_AB_RO2.csv".
%
%     NOTES:  1.  M-files decomp.m, parse_csv.m, rd_csv.m, and rd_roi6.m
%             must be in the current directory or path.
%
%             2.  See chk2_3d_2d.m for the program that just checks the
%             3D and 2D coordinates.
%
%     20-May-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Get Subject Directory and RO2 CSV File to Fix
%
cdirf = 'MRIR 07-BP\Visit2\T2S\Femur\';
csvf = '007_L_SAGAR_FEM_T2S_LD_V2_AB_RO2.csv';
%
% Read CSV File and Get Lines with Slices
%
irho = 5;
slsk = 33;
%
[flines,fhdrs,fids] = rd_csv(cdirf,csvf,slsk,irho);     % Femur
%
% Get Segmentation Points from Slice Lines
%
[fpts3d,fpts2d,~,fidn,fidx] = parse_csv(fhdrs,flines{fids});    % Femur
%
% Fix Coordinates by Removing "Bad" Point
%
npts = size(fpts3d,1);
idv = true(npts,1);
idv(19) = false;        % Remove 19th point
fpts3dro = fpts3d(idv,:);              % New corrected points
%
idv(19) = true;
idv(18) = false;
fpts2dro = fpts2d(idv,:);
%
% Use 3-D Coordinates for 2-D Points
%
fpts2dro = trnsf2pixel(fpts3dro,fpts2dro,fpts3dro);
%
% Update Segmentation File Lines with the New 3-D and 2-D Points
%
flines{fids} = line_upd(flines{fids},fidn,fidx,fpts3dro,fpts2dro);
%
% Write Updated Segmentation lines to a New CSV File 
%
write_csv(cdirf,csvf,flines);
%
return