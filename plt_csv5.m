%#######################################################################
%
%                     * PLoT CSV Files 5 Program *
%
%          M-File which reads the cartilage and bone digitization CSV
%     files and plots the data for visual verification.  The program
%     assumes the femur has three sets of digitizations (lateral,
%     medial and trochlea), the tibia has two sets of digitizations
%     (lateral and medial) and the patella has one set of digitizations.
%     Due to the angle of the knee alignment with the scanner, the
%     patella data is ignored.
%
%     NOTES:  1.  Matlab M-files plt_datsl.m and rd_roi6.m must be in
%             the current directory or path.
%
%             2.  Plots are output to PS files:  dig_plt5*.ps
%
%             3.  Femur digitization are assumed to be in subdirectory
%             FEMUR and tibia digitization in subdirectory TIBIA. 
%
%             4.  For the MRI reliability study subjects with left and
%             right and loaded and unloaded knees in each directory.
%             The CSV file names must contain the following keys:
%             Side:  _L_ for left knees or _R_ for right knees
%             Loading:  _LD_ for loaded knees or _UL_ for unloaded knees
%             Surface:  _SAG_ for subchondral bone or _SAGAR_ for
%                       articular cartilage
%             Bone:  _FEM_ for femur or _TIB_ for tibia (in separate
%                    subdirectories)
%
%             5.  CSV files with case insensitive "dup" in the file
%             names are ignored.
%
%     18-Jan-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Directories with Digitization CSV Files
%
dirstr = split(pwd,filesep);
dirstr = dirstr{end};
%
vdirs = ['Visit1'; 'Visit2'];          % Visit directories
nvdirs = size(vdirs,1);
%
mdirs = ['RHO'; 'T2S'];                % MRI type directories
nmdirs = size(mdirs,1);
%
% Side and Loading Codes in File Names
%
sides = ['*L*'; '*R*'];                % L - left and R - right
stxts = {'Left'; 'Right'};
%
loads = ['*LD*.csv'; '*UL*.csv'];      % LD - loaded and UL - unloaded
ltxts = {'Loaded'; 'Unloaded'};
%
% Loop Through Directories and Plot Digitizations
%
psfiles = 'dig_plt5';   % Initial PS file name
%
ido = false;            % Don't do patella
%
xnams = cell(64,1);     % List of CSV files
%
for k = 1:nvdirs
%
   vdir = vdirs(k,:);
   psv = vdir([1 6]);   % Visit number for PS file
   psfile = [psfiles '_' psv];
%
   for l = 1:nmdirs
%
      mdir = mdirs(l,:);
      psfil = [psfile '_' mdir '.ps'];
%
      for m = 1:2       % Side (L/R)
%
         side = sides(m,:);
%
         for n = 1:2    % Loads (Loaded/Unloaded)
%
            figure;
            orient landscape;
%
            ld = loads(n,:);
            idf = 32*k+16*l+8*m+4*n-56;
%
% Patella
% Needs updating for codes, left/right and loaded/unloaded knees
%
            if ido
              cdat = rd_roi6(fullfile(vdir,mdir,'*PAT_CART*.csv'));
              cdat = cdat.data;
              plt_datsl(cdat,'m.-');
              hold on;
              bdat = rd_roi6(fullfile(vdir,mdir,'*PAT_BONE*.csv'));
              bdat = bdat.data;
              plt_datsl(bdat,'k.-');
            end
%
% Femur
% 3 - TROCHLEA
%
            fdir = fullfile(vdir,mdir,'Femur');
            ffiles = fullfile(fdir,[side 'SAGAR_FEM' ld]);
            fnams = dir(ffiles);
            fnams = {fnams.name}';
            nfiles = size(fnams,1);
            if nfiles<1
              serr = sprintf([' *** ERROR in plt_csv5: No files', ...
                              'matching: %s'],ffiles);
              error(serr);
            elseif nfiles>1
              idx = ~contains(fnams,'dup','IgnoreCase',true);
              fnams = fnams(idx);
              nfiles = size(fnams,1);
              if nfiles<1
                serr = sprintf([' *** ERROR in plt_csv5: Not', ...
                                ' just one file matching: %s'],ffiles);
                error(serr);
              end
            end
            fnams = fullfile(fdir,fnams{1});
            xnams{idf-3} = fnams;
            fcdat = rd_roi6(fnams);
%
            nams = {fcdat.name}';
            fcdat1 = fcdat(1).data;
            plt_datsl(fcdat1,'r.-');
            fcdat2 = fcdat(2).data;
            plt_datsl(fcdat2,'b.-');
            fcdat3 = fcdat(3).data;
            plt_datsl(fcdat3,'m.-');
%
            ffiles = fullfile(fdir,[side 'SAG_FEM' ld]);
            fnams = dir(ffiles);
            fnams = {fnams.name}';
            nfiles = size(fnams,1);
            if nfiles<1
              serr = sprintf([' *** ERROR in plt_csv5: No files', ...
                              'matching: %s'],ffiles);
              error(serr);
            elseif nfiles>1
              idx = ~contains(fnams,'dup','IgnoreCase',true);
              fnams = fnams(idx);
              nfiles = size(fnams,1);
              if nfiles<1
                serr = sprintf([' *** ERROR in plt_csv5: Not', ...
                                ' just one file matching: %s'],ffiles);
                error(serr);
              end
            end
            fnams = fullfile(fdir,fnams{1});
            xnams{idf-2} = fnams;
            fbdat = rd_roi6(fnams);
%
            fbdat1 = fbdat(1).data;
            fbdat2 = fbdat(2).data;
            plt_datsl(fbdat1,'k.-');
            plt_datsl(fbdat2,'k.-');
            fbdat3 = fbdat(3).data;
            plt_datsl(fbdat3,'k.-');
%
% Tibia
%
            tdir = fullfile(vdir,mdir,'Tibia');
            tfiles = fullfile(tdir,[side 'SAGAR_TIB' ld]);
            fnams = dir(tfiles);
            fnams = {fnams.name}';
            nfiles = size(fnams,1);
            if nfiles<1
              serr = sprintf([' *** ERROR in plt_csv5: No files', ...
                              'matching: %s'],tfiles);
              error(serr);
            elseif nfiles>1
              idx = ~contains(fnams,'dup','IgnoreCase',true);
              fnams = fnams(idx);
              nfiles = size(fnams,1);
              if nfiles<1
                serr = sprintf([' *** ERROR in plt_csv5: Not', ...
                                ' just one file matching: %s'],tfiles);
                error(serr);
              end
            end
            fnams = fullfile(tdir,fnams{1});
            xnams{idf-1} = fnams;
            tcdat = rd_roi6(fnams);
%
            tcdat1 = tcdat(1).data;
            tcdat2 = tcdat(2).data;
            plt_datsl(tcdat1,'r.-');
            plt_datsl(tcdat2,'b.-');
%
            tfiles = fullfile(tdir,[side 'SAG_TIB' ld]);
            fnams = dir(tfiles);
            fnams = {fnams.name}';
            nfiles = size(fnams,1);
            if nfiles<1
              serr = sprintf([' *** ERROR in plt_csv5: No files', ...
                              'matching: %s'],tfiles);
              error(serr);
            elseif nfiles>1
              idx = ~contains(fnams,'dup','IgnoreCase',true);
              fnams = fnams(idx);
              nfiles = size(fnams,1);
              if nfiles<1
                serr = sprintf([' *** ERROR in plt_csv5: Not', ...
                                ' just one file matching: %s'],tfiles);
                error(serr);
              end
            end
            fnams = fullfile(tdir,fnams{1});
            xnams{idf} = fnams;
            tbdat = rd_roi6(fnams);
%
            tbdat1 = tbdat(1).data;
            tbdat2 = tbdat(2).data;
            plt_datsl(tbdat1,'k.-');
            plt_datsl(tbdat2,'k.-');
%
% Finish and Print Plot
%
            view(-110,15);
            axis equal;
            hx = xlabel('\leftarrow Lateral','FontSize',12, ...
                        'FontWeight','bold');
            hy = ylabel('Anterior \rightarrow','FontSize',12, ...
                        'FontWeight','bold');
            st = [dirstr ' ' vdir ' ' mdir ' ' stxts{m} ' ' ltxts{n}];
            title({st; ['Red-Lateral, Magenta-Groove and Blue-' ...
                   'Medial Cartilage']; ' Black-Bone'},'FontSize', ...
                   16,'FontWeight','bold','Interpreter','none');
            if m==1&&n==1
              print('-dpsc2','-r600','-fillpage',psfil);
            else
              print('-dpsc2','-r600','-fillpage','-append',psfil);
            end
%
%             pause;
%
         end            % n
%
      end               % m
%
      close all;
%
   end                  % l
%
end                     % k
%
xlsnam = 'dig_plt5.xlsx';
t = table(xnams,'VariableNames',{'Filenames'});
writetable(t,xlsnam);
%
return