%#######################################################################
%
%         * CONTACT CHecK of Segmentations MAY 2024 Program *
%
%          M-File which reads the segmentation CSV files to check if
%     the femur and tibia segmentations overlap.  If they do overlap,
%     the program measures the maximum amount of overlap.
%
%          The slices with overlap are plotted to Postscript files
%     contact_chk*.ps in the Results\ContactChecksMay24 folder.  The
%     maximum overlaps are also saved to the MS-Excel spreadsheet, 
%     contact_chkMay24.xlsx, in the Results\ContactChecksMay24 folder.
%
%     NOTES:  1.  M-files lsect3.m, lsect4.m, lsect5.m, mxd2lins.m,
%             pts2lin.m, rd_roi6.m and rd_roiso.m must be in the
%             current directory or path.
%
%     15-May-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results','ContactChecksMay24');    % Results directory
%
ifirst = true;          % First write to file
xlsnam = 'contact_chkMay24.xlsx';      % Results spreadsheet
xlsnam = fullfile(resdir,xlsnam);      % Include output directory
hdrs = {'Subject' 'Visit' 'Result' 'Leg' 'Load' 'Comprt' 'Slice', ...
        'MaxDist'};
%
psnam = fullfile(resdir,'contact_chkMay24_');    % Start of PS file name
pstyp = '.ps';          % PS file type
%
% Get Subject Directories and Visit Subdirectories
%
sdirs = dir('MRIR *');
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
nsubj = size(sdirs,1);
%
vdirs = {'Visit1'; 'Visit2'};          % Visit directories
nvisit = size(vdirs,1);
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
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit
%   Index 3 - Result - 1 = T1rho and 2 = T2*
%   Index 4 - Leg - 1 = left and 2 = right
%   Index 5 - Load - 1 = unloaded and 2 = loaded
%   Index 6 - Compartment - 1 = lateral and 2 = medial
%
dmx = cell(nsubj,nvisit,2,2,2,2);      % Maximum distances
sls = cell(nsubj,nvisit,2,2,2,2);      % Slice numbers
dmxs = zeros(nsubj,nvisit,2,2,2,2);    % Maximum distances for leg, load and compartment
%
% T1rho Segmentations
%
rdir = 'RHO';
ires = 1;               % T1rho
itroch = true;
%
% Loop through Subjects
%
for ks = 1:nsubj
%
% Get Subject Directory, Name and Number
%
   sdir = sdirs{ks};                   % Current subject directory
   subjnam = sdir(6:end);              % Subject name as text
   subj = eval(subjnam(1:2));          % Subject number
%
   psnams = [psnam subjnam];           % Add subject to PS file name
%
% Loop through Visits
%
   for kv = 1:nvisit
%
% Get Visit Subdirectory, Name and Number
%
      vdir = vdirs{kv};                % Current visit directory
      vstr = int2str(kv);              % Visit number as text
      vnam = ['Visit ' vstr];          % Visit name as text
      vid = kv-1;                      % Visit number
%
      psnamv = [psnams '_V' vstr];     % Add visit to PS file name
%
% T1rho Identifier
%
      psnamr = [psnamv '_T1R_'];       % Add result type to PS file name
%
% Directory with Data Matlab MAT Files
%
      rdirk = fullfile(sdir,vdir,rdir);     % Directory with data
%
% Loop through Legs
%
      for kl = 1:2
%
% Get Leg
%
         leg = legs(kl);
         lnam = lnams{kl};
%
% Loop through Loads
%
         for l = 1:2
%
% Get Load
%
            ld = lds(l,:);
            ldnam = ldnams{l};
%
% Add Leg and Load to PS File Name
%
            psnamf = [psnamr leg '_' ld pstyp];  % Add leg and load to PS file name
            iplt1 = true;              % First plot in new PS file
%
% Read Segmentations
%
            brois = rd_roiso(rdirk,leg,ld,itroch,4);
%
% Get Femur Data
%
            fc = brois(1).rois(1);     % Femur cartilage ROIs
            fsl = fc.slice;            % All femur slices
            fdat = [fc.roi.data]';     % All femur data
%
% Tibia Compartment Data
%
            for m = 1:2 % 1 - lateral and 2 - medial
%
               tc = brois(2).rois(1);  % Tibia cartilage ROIs
%
               tcc = tc.roi(m);        % Tibia compartment
%                tcl = tc.roi(1);        % Tibia lateral compartment
%                tcm = tc.roi(2);        % Tibia medial compartment
               tsl = tcc.imageno;      % Compartment slices
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
                       idcs(mi-1:mi,:) = sort(idcs(mi-1:mi,:)); % Make sure both lines are going in the same direction
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
                    end
%
                    if iplt1
                      print('-dpsc2','-r600','-fillpage',psnamf);
                      iplt1 = false;
                    else
                      print('-dpsc2','-r600','-fillpage','-append', ...
                            psnamf);
                    end
                    close;
%
% Save Slice Maximum Overlap and Slice Number
%
                    dmx{ks,kv,ires,kl,l,m} = [dmx{ks,kv,ires,kl,l, ...
                                              m}; dd];
                    sls{ks,kv,ires,kl,l,m} = [sls{ks,kv,ires,kl,l, ...
                                              m}; fsl(idf(n))];
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
               end      % Slices loop - n
%
               if ~isempty(dmx{ks,kv,ires,kl,l,m})
                 dmxs(ks,kv,ires,kl,l,m) = max(dmx{ks,kv,ires,kl,l,m});
               end
%
            end         % End of tibia compartment loop - m
%
         end            % End of load loop - l
%
      end               % End of leg loop - kl
%
   end                  % End of visit loop - kv
%
end                     % End of subject loop - ks
%
% T2* Segmentations
%
rdir = 'T2S';
ires = 2;               % T2*
%
% Loop through Subjects
%
for ks = 1:nsubj
%
% Get Subject Directory, Name and Number
%
   sdir = sdirs{ks};                   % Current subject directory
   subjnam = sdir(6:end);              % Subject name as text
   subj = eval(subjnam(1:2));          % Subject number
%
   psnams = [psnam subjnam];           % Add subject to PS file name
%
% Loop through Visits
%
   for kv = 1:nvisit
%
% Get Visit Subdirectory, Name and Number
%
      vdir = vdirs{kv};                % Current visit directory
      vstr = int2str(kv);              % Visit number as text
      vnam = ['Visit ' vstr];          % Visit name as text
      vid = kv-1;                      % Visit number
%
      psnamv = [psnams '_V' vstr];     % Add visit to PS file name
%
% T2* Identifier
%
      psnamr = [psnamv '_T2S_'];       % Add result type to PS file name
%
% Directory with Data Matlab MAT Files
%
      rdirk = fullfile(sdir,vdir,rdir);     % Directory with data
%
% Loop through Legs
%
      for kl = 1:2
%
% Get Leg
%
         leg = legs(kl);
         lnam = lnams{kl};
%
% Loop through Loads
%
         for l = 1:2
%
% Get Load
%
            ld = lds(l,:);
            ldnam = ldnams{l};
%
% Add Leg and Load to PS File Name
%
            psnamf = [psnamr leg '_' ld pstyp];  % Add leg and load to PS file name
            iplt1 = true;              % First plot in new PS file
%
% Read Segmentations
%
            brois = rd_roiso(rdirk,leg,ld,itroch,5);
%
% Get Femur Data
%
            fc = brois(1).rois(1);     % Femur cartilage ROIs
            fsl = fc.slice;            % All femur slices
            fdat = [fc.roi.data]';     % All femur data
%
% Tibia Compartment Data
%
            for m = 1:2 % 1 - lateral and 2 - medial
%
               tc = brois(2).rois(1);  % Tibia cartilage ROIs
%
               tcc = tc.roi(m);        % Tibia compartment
%                tcl = tc.roi(1);        % Tibia lateral compartment
%                tcm = tc.roi(2);        % Tibia medial compartment
               tsl = tcc.imageno;      % Compartment slices
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
                    title({['Subject ' subjnam ' ' vnam ' T2*']; ...
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
                       idcs(mi-1:mi,:) = sort(idcs(mi-1:mi,:)); % Make sure both lines are going in the same direction
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
                    end
%
                    if iplt1
                      print('-dpsc2','-r600','-fillpage',psnamf);
                      iplt1 = false;
                    else
                      print('-dpsc2','-r600','-fillpage','-append', ...
                            psnamf);
                    end
                    close;
%
% Save Slice Maximum Overlap and Slice Number
%
                    dmx{ks,kv,ires,kl,l,m} = [dmx{ks,kv,ires,kl,l, ...
                                              m}; dd];
                    sls{ks,kv,ires,kl,l,m} = [sls{ks,kv,ires,kl,l, ...
                                              m}; fsl(idf(n))];
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
               end      % Slices loop - n
%
               if ~isempty(dmx{ks,kv,ires,kl,l,m})
                 dmxs(ks,kv,ires,kl,l,m) = max(dmx{ks,kv,ires,kl,l,m});
               end
%
            end         % End of tibia compartment loop - m
%
         end            % End of load loop - l
%
      end               % End of leg loop - kl
%
   end                  % End of visit loop - kv
%
end                     % End of subject loop - ks
%
% Get Overlaps by Analysis Type
%
dmx1 = squeeze(dmx(:,:,1,:,:,:));      % T1rho
dmx2 = squeeze(dmx(:,:,2,:,:,:));      % T2*
%
dmx1 = cell2mat({dmx1{:}}');
dmx2 = cell2mat({dmx2{:}}');
%
dmx1 = dmx1(dmx1>0);
dmx2 = dmx2(dmx2>0);
%
% Overlap Descriptive Statistics
%
n1 = size(dmx1,1);
n2 = size(dmx2,1);
%
n21 = nnz(dmx1>2);      % Overlaps > 2 pixels
n22 = nnz(dmx2>2);      % Overlaps > 2 pixels
%
p21 = 100*n21/n1;
p22 = 100*n22/n2;
%
mn1 = mean(dmx1);
mx1 = max(dmx1);
sd1 = std(dmx1);
mn2 = mean(dmx2);
mx2 = max(dmx2);
sd2 = std(dmx2);
%
% Overlap Histograms
%
figure;
orient landscape;
histogram(dmx1);
xlabel('Overlap (pixels)','FontSize',12,'FontWeight','bold');
ylabel('Frequency','FontSize',12,'FontWeight','bold');
title('T1\rho Overlaps','FontSize',16,'FontWeight','bold');
%
axlim = axis;
hold on;
plot([2 2],axlim(3:4),'k:','LineWidth',0.5);
ptxt{3} = sprintf('SD = %.2f pixels',sd1);
ptxt{2} = sprintf('Maximum = %.2f pixels',mx1);
ptxt{1} = sprintf('Mean = %.2f pixels',mn1);
text(0.9,sum(axlim(3:4))/2,ptxt,'FontSize',12,'FontWeight','bold');
ptxt = sprintf('%.2f%% of overlaps > 2 pixels',p21);
ptxt = {[int2str(n21) ' overlaps > 2 pixels']; ptxt};
text(2.1,sum(axlim(3:4))/4,ptxt,'FontSize',12,'FontWeight','bold');
%
psnamh = fullfile(resdir,'overlap_hist.ps');
print('-dpsc2','-r600','-fillpage',psnamh);
%
figure;
orient landscape;
histogram(dmx2);
xlabel('Overlap (pixels)','FontSize',12,'FontWeight','bold');
ylabel('Frequency','FontSize',12,'FontWeight','bold');
title('T2* Overlaps','FontSize',16,'FontWeight','bold');
%
axlim = axis;
hold on;
plot([2 2],axlim(3:4),'k:','LineWidth',0.5);
ptxt{3} = sprintf('SD = %.2f pixels',sd2);
ptxt{2} = sprintf('Maximum = %.2f pixels',mx2);
ptxt{1} = sprintf('Mean = %.2f pixels',mn2);
text(0.9,sum(axlim(3:4))/2,ptxt,'FontSize',12,'FontWeight','bold');
ptxt = sprintf('%.2f%% of overlaps > 2 pixels',p22);
ptxt = {[int2str(n22) ' overlaps > 2 pixels']; ptxt};
text(2.1,sum(axlim(3:4))/4,ptxt,'FontSize',12,'FontWeight','bold');
print('-dpsc2','-r600','-append','-fillpage',psnamh);
%
return