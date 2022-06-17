%#######################################################################
%
%             * CONTACT CHecK of Segmentations 3 Program *
%
%          M-File which reads the modified segmentation CSV files from
%     a list in a MS-Excel spreadsheet and checks if the femur and
%     tibia segmentations overlap.  If they do overlap, the program
%     measures the maximum amount of overlap.
%
%          The slices with overlap are plotted to Postscript files
%     contact_chk3_*.ps in the Results\ContactChecks2 folder.  The
%     maximum overlaps are also saved to the MS-Excel spreadsheet, 
%     contact_chk3.xlsx, in the Results\ContactChecks2 folder.
%
%     NOTES:  1.  M-files lsect3.m, lsect4.m, lsect5.m, mxd2lins.m,
%             pts2lin.m, rd_roi6.m, rd_rois.m and rd_rois2.m must be in
%             the current directory or path.
%
%             2.  MS-Excel spreadsheet "Modified_ROIs_09Jun2022.xlsx"
%             must be in the subdirectory "Results\NACOB_Final\".
%
%     14-Jun-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results','ContactChecks3');   % Results directory
%
ifirst = true;          % First write to file
% xlsnam = 'contact_chk.xlsx';           % Results spreadsheet
% xlsnam = 'contact_chk_RO.xlsx';        % Results spreadsheet
% xlsnam = 'contact_chk2.xlsx';          % Results spreadsheet
xlsnam = 'contact_chk3.xlsx';          % Results spreadsheet
xlsnam = fullfile(resdir,xlsnam);      % Include output directory
hdrs = {'Subject' 'Visit' 'Result' 'Leg' 'Load' 'Comprt' 'Slice', ...
        'MaxDist'};
%
psnam = fullfile(resdir,'contact_chk3_');   % Start of PS file name
pstyp = '.ps';          % PS file type
%
% Get Modified CSV File Names
%
rdir = 'Results\NACOB_Final';
xlsnam_inp = 'Modified_ROIs_09Jun2022.xlsx';
t = readtable(fullfile(rdir,xlsnam_inp));
idx = t.NeedFix;
idx = find(idx==1);
csvfiles = t.CSVFiles;
csvfiles = csvfiles(idx);
nfiles = size(csvfiles,1);
%
% Get Subject and Visit Directories
%
sdirs = dir('MRIR*');
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
%
vdirs = {'Visit1'; 'Visit2'};          % Visit directories
%
% Get Legs and Loads Variables
%
legs = ['L'; 'R'];      % Left/Right
lnams = {'Left Leg'; 'Right Leg'};
%
lds = ['UL'; 'LD'];     % Unloaded/loaded
ldnams = {'Unloaded'; 'Loaded'};
%
cmprt  = {'Lateral'; 'Medial'};
%
% Initialize Results Variables
%
dmx = cell(nfiles,2);   % Maximum distances
sls = cell(nfiles,2);   % Slice numbers
dmxs = zeros(nfiles,2); % Maximum distances for leg, load and compartment
%
itroch = true;          % Read trochlear slices
%
% Loop through CSV Files
%
for kf = 1:nfiles
%
% Get CSV File
%
   csvfile = csvfiles{kf};
%
% Get Subject Directory, Name and Number
%
   ids = contains(sdirs,csvfile(1:2));
   sdir = sdirs{ids};                  % Current subject directory
   subjnam = sdir(6:end);              % Subject name as text
   subj = eval(subjnam(1:2));          % Subject number
%
   psnams = [psnam subjnam];           % Add subject to PS file name
%
% Get Visit Subdirectory, Name and Number
%
   vdir = csvfile(4:9);                % Current visit directory
   vstr = csvfile(9);                  % Visit number as text
   vnam = ['Visit ' vstr];             % Visit name as text
   vid = eval(csvfile(9))-1;           % Visit number
%
   psnamv = [psnams '_V' vstr];        % Add visit to PS file name
%
% Analysis Segmentations
%
   if contains(csvfile,'RHO')
     rdir = 'RHO';
     ires = 1;          % T1rho
     irho = 4;
   else
     rdir = 'T2S';
     ires = 2;          % T2*
     irho = 5;
   end
%
   psnamr = [psnamv '_' rdir '_'];     % Add result type to PS file name
%
% Directory with CSV Segmentations Files
%
   rdirk = fullfile(sdir,vdir,rdir);   % Directory with data
%
% Get Leg
%
   kl = contains(csvfile,'_R_')+1;
%
% Get Leg
%
   leg = legs(kl);
   lnam = lnams{kl};
%
% Get Loads
%
   l = contains(csvfile,'_LD_')+1;
%
% Get Load
%
   ld = lds(l,:);
   ldnam = ldnams{l};
%
% Add Leg and Load to PS File Name
%
   psnamf = [psnamr leg '_' ld pstyp]; % Add leg and load to PS file name
   iplt1 = true;             % First plot in new PS file
%
% Read Segmentations
%
   brois = rd_rois2(rdirk,leg,ld,itroch,irho);
%
% Get Femur Data
%
   fc = brois(1).rois(1);              % Femur cartilage ROIs
   fsl = fc.slice;                     % All femur slices
   fdat = [fc.roi.data]';              % All femur data
%
% Tibia Compartment Data
%
   for m = 1:2          % 1 - lateral and 2 - medial
%
      tc = brois(2).rois(1);           % Tibia cartilage ROIs
%
      tcc = tc.roi(m);                 % Tibia compartment
      tsl = tcc.imageno;               % Compartment slices
%
% Find Slices with Both Tibia and Femur Segmentations
%
      [~,idt,idf] = intersect(tsl,fsl);
%
% Loop through Compartment Slices
%
      ns = size(idt,1);
%
      for n = 1:ns
%
         fdats = fdat{idf(n)};
         tdats = tcc.data{idt(n)};
%
% Find Any Intersections Between the Tibia and Femur Segmentations
%
         [ipts,~,idcs] = lsect5(tdats,fdats);
%
         if ~isempty(ipts)
           figure;
           orient landscape;
%
           plot(fdats(:,1),fdats(:,2),'b.-','LineWidth', ...
                1.5,'MarkerSize',12);
           hold on;
           plot(tdats(:,1),tdats(:,2),'g.-','LineWidth', ...
                1.5,'MarkerSize',12);
           plot(ipts(:,1),ipts(:,2),'ro','LineWidth',1.5, ...
                'MarkerSize',8);
           axis equal;
           set(gca,'XDir','reverse');
           set(gca,'YDir','reverse');
%
           title({['Subject ' subjnam ' ' vnam ' T1\rho']; ...
                  [cmprt{m} ' ' lnam ', ' ldnam ', Slice ' ...
                  int2str(tsl(idt(n)))]},'FontSize',16, ...
                 'FontWeight','bold');
%
           ni = size(ipts,1);
           dd = 0;
%
           for km = 1:ni/2
              mi = 2*km;
              idfx = fdats(:,1)>ipts(mi-1,1)&fdats(:,1)< ...
                           ipts(mi,1);
              idtx = tdats(:,1)>ipts(mi-1,1)&tdats(:,1)< ...
                           ipts(mi,1);
              nf = nnz(idfx);
              nt = nnz(idtx);
              df = 0;
              dt = 0;
%
              idcs(mi-1:mi,:) = sort(idcs(mi-1:mi,:));     % Make sure both lines are going in the same direction
%
              if nf>0
                idtl = idcs(mi-1,1):idcs(mi,1)+1;
                txyz = tdats(idtl,:);
                fxyz = fdats(idfx,:);
                [df,idfmx,xyzf] = mxd2lins(txyz,fxyz);
                if idfmx>0
                  plot([fxyz(idfmx,1); xyzf(:,1)], ...
                       [fxyz(idfmx,2); xyzf(:,2)],'r-', ...
                       'LineWidth',2);
                  text(xyzf(:,1),xyzf(:,2),sprintf('%.2f',df));
                end
              end
%
              if nt>0
                idfl = idcs(mi-1,2):idcs(mi,2)+1;
                txyz = tdats(idtx,:);
                fxyz = fdats(idfl,:);
                [dt,idtmx,xyzt] = mxd2lins(fxyz,txyz);
                if idtmx>0
                  plot([txyz(idtmx,1); xyzt(:,1)], ...
                       [txyz(idtmx,2); xyzt(:,2)],'m-', ...
                       'LineWidth',2);
                  text(xyzt(:,1),xyzt(:,2),sprintf('%.2f',dt));
                end
              end
%
% Get Slice Maximum Overlap
%
              if df>dt&&df>dd
                dd = df;
              elseif dt>dd
                dd = dt;
              end
%
            end         % if ~isempty(ipts)
%
            if iplt1
              print('-dpsc2','-r600','-fillpage',psnamf);
              iplt1 = false;
            else
              print('-dpsc2','-r600','-fillpage','-append',psnamf);
            end
            close;
%
% Save Slice Maximum Overlap and Slice Number
%
            dmx{kf,m} = [dmx{kf,m}; dd];
            sls{kf,m} = [sls{kf,m}; fsl(idf(n))];
%
% Create and Write Table of Results
%
            dat = [subj vid ires-1 kl-1 l-1 m-1 fsl(idf(n)) dd];
            t = array2table(dat,'VariableNames',hdrs);
%
            if ifirst
              writetable(t,xlsnam,'WriteMode','replacefile');
              ifirst = false;
            else
              writetable(t,xlsnam,'WriteMode','append', ...
                         'WriteVariableNames',false);
            end
%       
         end   % If femur and tibia segmentation intersect?
%
      end               % Slices loop - n
%
      if ~isempty(dmx{kf,m})
        dmxs(kf,m) = max(dmx{kf,m});
      end
%
   end                  % End of tibia compartment loop - m
%
end                     % End of CSV file loop - kf
%
return