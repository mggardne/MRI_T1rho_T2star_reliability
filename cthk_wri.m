%#######################################################################
%
%                * Cartilage THicKness WRIte Program *
%
%          M-File which reads the cartilage thicknesses for the MRI
%     reliability study femurs and tibias from the MS-Excel
%     spreadsheets "Femoral Thicknesses.xlsx" and 
%     "Tibial Thicknesses.xlsx".  The femoral cartilage thicknesses are
%     written to six MS-Excel spreadsheets, fem_lat_ant.xlsx,
%     fem_lat_ctr.xlsx, fem_lat_pos.xlsx, fem_med_ant.xlsx,
%     fem_med_ctr.xlsx, and fem_med_pos.xlsx.  The tibial cartilage
%     thicknesses are written to six MS-Excel spreadsheets,
%     tib_lat_ant.xlsx, tib_lat_ctr.xlsx, tib_lat_pos.xlsx,
%     tib_med_ant.xlsx, tib_med_ctr.xlsx, and tib_med_pos.xlsx.
%
%     NOTE:   1.  The MS-Excel spreadsheets "Femoral Thicknesses.xlsx"
%             and "Tibial Thicknesses.xlsx" must be in the current
%             directory.
%
%     07-Aug-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Regions of Interest (ROIs) and Grid Data
%
nr = 6;                 % Number of regions
nc = 2;                 % Number of compartments
nd = 3;                 % Number of divisions
%
cnams = ['lat'; 'med'];                % Compartment names
dnams = ['pos'; 'ctr'; 'ant'];         % Division names
%
nt = 3;                 % Number of MRI scan types (FFE, T1rho, and T2*)
%
p = 'Pt_';              % Prefix for grid point numbers
%
% Get Femur Cartilage Thicknesses
%
xlsnami = 'Femoral Thicknesses.xlsx';  % Input femur spreadsheet filename
[~,shts] = xlsfinfo(xlsnami);
id0 = startsWith(shts,'0');
shts = shts(id0)';      % Name of subject sheets in the spreadsheet
%
nsubj = size(shts,1);   % Number of subjects (10)
subjn = char(shts);
subjn = str2double(cellstr(subjn(:,1:3)));       % Subject numbers
%
% Loop through Subjects (Sheets)
%
xlsnamos = 'fem_';      % Initial output femur spreadsheet filenames
xlsfext = '.xlsx';      % Output spreadsheet filename extension
%
for ks = 1:nsubj        % Loop through subjects (sheets)
%
   opts = detectImportOptions(xlsnami,'Sheet',shts{ks});
   ncol = size(opts.VariableTypes,2);
   opts.VariableTypes = cellstr(repmat('double',ncol,1))'; % Read numbers
%
   t = readtable(xlsnami,opts);        % Read spreadsheet sheet
%
% Loop through Table Columns
%
   for kt = 1:nt        % Loop through MRI scan types
%
      for kd = 1:nd     % Loop through divisions
%
         for kc = 1:nc  % Loop through compartments
%
% Read Cartilage Thicknesses and Grid Point Numbers
%
            idr = kd*nc+kc-nc;         % Index to region
            icol = 13*kt+2*idr-14;     % Index to table column
%
            idv = ~isnan(t{:,icol+1}); % Get valid grid points
            grdpts = t{idv,icol+1};    % Get grid points
            thks = t{idv,icol};        % Get thicknesses
%
% Write Cartilage Thicknesses and Grid Point Numbers
%
            npts = sum(idv);           % Number of valid grid points
            pt_lbls = mat2cell([repmat('Pt_',npts,1) ...
                                int2str((grdpts))],ones(npts,1))';
            pt_lbls = strrep(pt_lbls,' ','');    % Grid point labels
            lbls = ['Subject' 'ScanType' 'Compartment' 'Division' ...
                    pt_lbls];          % Labels for analysis variables
%
            id = [subjn(ks) kt-1 kc-1 3-kd];
            tout = array2table([id thks'],'VariableNames',lbls);
%
            xlsnamo = [xlsnamos cnams(kc,:) '_' dnams(kd,:) xlsfext];
%
            if ks==1&&kt==1
              writetable(tout,xlsnamo,'WriteMode','replacefile');

            else
              writetable(tout,xlsnamo,'WriteMode','append', ...
                         'WriteVariableNames',false);
            end
%
         end
%
      end
%
   end
%
end
%
% Get Tibia Cartilage Thicknesses
%
% Note:  Tibia sheets and subject numbers are the same as the femur.
%
xlsnami = 'Tibial Thicknesses.xlsx';   % Input tibial spreadsheet filename
%
% Loop through Subjects (Sheets)
%
xlsnamos = 'tib_';      % Initial output tibia spreadsheet filenames
xlsfext = '.xlsx';      % Output spreadsheet filename extension
%
for ks = 1:nsubj        % Loop through subjects (sheets)
%
   opts = detectImportOptions(xlsnami,'Sheet',shts{ks});
   ncol = size(opts.VariableTypes,2);
   opts.VariableTypes = cellstr(repmat('double',ncol,1))'; % Read numbers
%
   t = readtable(xlsnami,opts);        % Read spreadsheet sheet
%
% Loop through Table Columns
%
   for kt = 1:nt        % Loop through MRI scan types
%
      for kd = 1:nd     % Loop through divisions
%
         for kc = 1:nc  % Loop through compartments
%
% Read Cartilage Thicknesses and Grid Point Numbers
%
            idr = kd*nc+kc-nc;         % Index to region
            icol = 13*kt+2*idr-14;     % Index to table column
%
            idv = ~isnan(t{:,icol+1}); % Get valid grid points
            grdpts = t{idv,icol+1};    % Get grid points
            thks = t{idv,icol};        % Get thicknesses
%
% Write Cartilage Thicknesses and Grid Point Numbers
%
            npts = sum(idv);           % Number of valid grid points
            pt_lbls = mat2cell([repmat('Pt_',npts,1) ...
                                int2str((grdpts))],ones(npts,1))';
            pt_lbls = strrep(pt_lbls,' ','');    % Grid point labels
            lbls = ['Subject' 'ScanType' 'Compartment' 'Division' ...
                    pt_lbls];          % Labels for analysis variables
%
            id = [subjn(ks) kt-1 kc-1 3-kd];
            tout = array2table([id thks'],'VariableNames',lbls);
%
            xlsnamo = [xlsnamos cnams(kc,:) '_' dnams(kd,:) xlsfext];
%
            if ks==1&&kt==1
              writetable(tout,xlsnamo,'WriteMode','replacefile');

            else
              writetable(tout,xlsnamo,'WriteMode','append', ...
                         'WriteVariableNames',false);
            end
%
         end
%
      end
%
   end
%
end
%
return