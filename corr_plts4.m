%
% Subject Numbers or Dots for Markers?
%
imrk = true;            % Use subject numbers
% imrk = false;           % Use dots
%
% Results Directory
%
rdir = 'Results\NACOB_Final4';
%
% T1rho Means
%
% mri_fitr4stat.xlsx has results from poor registrations set to -1
%
num = xlsread(fullfile(rdir,'mri_fitr4stat.xlsx'),'T1rho');
%
% Indices
%
ivis = logical(num(:,2));              % Visit
ileg = logical(num(:,4));              % Leg
ild = logical(num(:,5));               % Load
icmprt = logical(num(:,6));            % Compartment
ibone = logical(num(:,7));             % Bone
t1r = num(:,12);                       % T1rho
%
% Change Invalid T1rho from -1 to NaN
%
t1r(t1r==-1) = NaN;
%
% Get Subject Numbers
%
subj = unique(num(:,1));
subjt = cellstr(int2str(subj));
subjt = strrep(subjt,' ','');
%
% Get Indices for Visits 1 and 2
%
ilegv(:,2) = ileg(ivis);
ilegv(:,1) = ileg(~ivis);
ildv(:,2) = ild(ivis);
ildv(:,1) = ild(~ivis);
icmprtv(:,2) = icmprt(ivis);
icmprtv(:,1) = icmprt(~ivis);
ibonev(:,2) = ibone(ivis);
ibonev(:,1) = ibone(~ivis);
t1rv(:,2) = t1r(ivis);
t1rv(:,1) = t1r(~ivis);
%
% Get Indices for Left and Right Legs
%
% Indices:
% 1 - Visit 1 Left Leg
% 2 - Visit 1 Right Leg
% 3 - Visit 2 Left Leg
% 4 - Visit 2 Right Leg
%
ildvl(:,4) = ildv(ilegv(:,2),2);
ildvl(:,3) = ildv(~ilegv(:,2),2);
ildvl(:,2) = ildv(ilegv(:,1),1);
ildvl(:,1) = ildv(~ilegv(:,1),1);
%
icmprtvl(:,4) = icmprtv(ilegv(:,2),2);
icmprtvl(:,3) = icmprtv(~ilegv(:,2),2);
icmprtvl(:,2) = icmprtv(ilegv(:,1),1);
icmprtvl(:,1) = icmprtv(~ilegv(:,1),1);
%
ibonevl(:,4) = ibonev(ilegv(:,2),2);
ibonevl(:,3) = ibonev(~ilegv(:,2),2);
ibonevl(:,2) = ibonev(ilegv(:,1),1);
ibonevl(:,1) = ibonev(~ilegv(:,1),1);
%
t1rvl(:,4) = t1rv(ilegv(:,2),2);
t1rvl(:,3) = t1rv(~ilegv(:,2),2);
t1rvl(:,2) = t1rv(ilegv(:,1),1);
t1rvl(:,1) = t1rv(~ilegv(:,1),1);
%
% Get Indices for Unloaded and Loaded
%
% Indices:
% 1 - Visit 1 Left Leg Unloaded
% 2 - Visit 1 Left Leg Loaded
% 3 - Visit 1 Right Leg Unloaded
% 4 - Visit 1 Right Leg Loaded
% 5 - Visit 2 Left Leg Unloaded
% 6 - Visit 2 Left Leg Loaded
% 7 - Visit 2 Right Leg Unloaded
% 8 - Visit 2 Right Leg Loaded
%
for k = 8:-2:2
   l = k/2;
   icmprtvll(:,k) = icmprtvl(ildvl(:,l),l);
   ibonevll(:,k) = ibonevl(ildvl(:,l),l);
   t1rvll(:,k) = t1rvl(ildvl(:,l),l);
end
%
for k = 7:-2:1
   l = (k+1)/2;
   icmprtvll(:,k) = icmprtvl(~ildvl(:,l),l);
   ibonevll(:,k) = ibonevl(~ildvl(:,l),l);
   t1rvll(:,k) = t1rvl(~ildvl(:,l),l);
end
%
% Get Indices for Lateral and Medial Compartments
%
% Indices:
% 1 - Visit 1 Left Leg Unloaded Lateral
% 2 - Visit 1 Left Leg Unloaded Medial
% 3 - Visit 1 Left Leg Loaded Lateral
% 4 - Visit 1 Left Leg Loaded Medial
% 5 - Visit 1 Right Leg Unloaded Lateral
% 6 - Visit 1 Right Leg Unloaded Medial
% 7 - Visit 1 Right Leg Loaded Lateral
% 8 - Visit 1 Right Leg Loaded Medial
% 9 - Visit 2 Left Leg Unloaded Lateral
% 10 - Visit 2 Left Leg Unloaded Medial
% 11 - Visit 2 Left Leg Loaded Lateral
% 12 - Visit 2 Left Leg Loaded Medial
% 13 - Visit 2 Right Leg Unloaded Lateral
% 14 - Visit 2 Right Leg Unloaded Medial
% 15 - Visit 2 Right Leg Loaded Lateral
% 16 - Visit 2 Right Leg Loaded Medial
%
for k = 16:-2:2
   l = k/2;
   ibonevllc(:,k) = ibonevll(icmprtvll(:,l),l);
   t1rvllc(:,k) = t1rvll(icmprtvll(:,l),l);
end
%
for k = 15:-2:1
   l = (k+1)/2;
   ibonevllc(:,k) = ibonevll(~icmprtvll(:,l),l);
   t1rvllc(:,k) = t1rvll(~icmprtvll(:,l),l);
end
%
% Get Indices for the Femur and Tibia
%
% Indices:
% 1 - Visit 1 Left Leg Unloaded Lateral Femur
% 2 - Visit 1 Left Leg Unloaded Lateral Tibia
% 3 - Visit 1 Left Leg Unloaded Medial Femur
% 4 - Visit 1 Left Leg Unloaded Medial Tibia
% 5 - Visit 1 Left Leg Loaded Lateral Femur
% 6 - Visit 1 Left Leg Loaded Lateral Tibia
% 7 - Visit 1 Left Leg Loaded Medial Femur
% 8 - Visit 1 Left Leg Loaded Medial Tibia
% 9 - Visit 1 Right Leg Unloaded Lateral Femur
% 10 - Visit 1 Right Leg Unloaded Lateral Tibia
% 11 - Visit 1 Right Leg Unloaded Medial Femur
% 12 - Visit 1 Right Leg Unloaded Medial Tibia
% 13 - Visit 1 Right Leg Loaded Lateral Femur
% 14 - Visit 1 Right Leg Loaded Lateral Tibia
% 15 - Visit 1 Right Leg Loaded Medial Femur
% 16 - Visit 1 Right Leg Loaded Medial Tibia
% 17 - Visit 2 Left Leg Unloaded Lateral Femur
% 18 - Visit 2 Left Leg Unloaded Lateral Tibia
% 19 - Visit 2 Left Leg Unloaded Medial Femur
% 20 - Visit 2 Left Leg Unloaded Medial Tibia
% 21 - Visit 2 Left Leg Loaded Lateral Femur
% 22 - Visit 2 Left Leg Loaded Lateral Tibia
% 23 - Visit 2 Left Leg Loaded Medial Femur
% 24 - Visit 2 Left Leg Loaded Medial Tibia
% 25 - Visit 2 Right Leg Unloaded Lateral Femur
% 26 - Visit 2 Right Leg Unloaded Lateral Tibia
% 27 - Visit 2 Right Leg Unloaded Medial Femur
% 28 - Visit 2 Right Leg Unloaded Medial Tibia
% 29 - Visit 2 Right Leg Loaded Lateral Femur
% 30 - Visit 2 Right Leg Loaded Lateral Tibia
% 31 - Visit 2 Right Leg Loaded Medial Femur
% 32 - Visit 2 Right Leg Loaded Medial Tibia
%
for k = 32:-2:2
   l = k/2;
   t1rvllcb(:,k) = t1rvllc(ibonevllc(:,l),l);
end
%
for k = 31:-2:1
   l = (k+1)/2;
   t1rvllcb(:,k) = t1rvllc(~ibonevllc(:,l),l);
end
%
% Labels and Colors
%
lleg = {'Left'; 'Right'};
lld = {'Unloaded'; 'Loaded'};
lcmprt = {'Lateral'; 'Medial'};
lbone = {'Femur'; 'Tibia'};
%
colr = [0 0.2 0.8; 0 0.5 0];               % Blue and Green
%
% Plot Leg and Load
%
mnmx = [min(min(t1rvll)) max(max(t1rvll))];
%
hl = gobjects(2,1);     % Line handles
%
for kl = 1:2            % Leg
%
   l = 2*kl;
   figure;
   orient landscape;
   hold on;
%
   for kd = 1:2         % Load
%
      l1 = l+kd-2;
      l2 = l1+4;
%
      hl(kd) = plot(t1rvll(:,l1),t1rvll(:,l2),'.','Color', ...
                    colr(kd,:),'MarkerSize',8,'LineWidth',1);
      idv = ~isnan(t1rvll(:,l1));      % NaNs only on Visit 1
      [b,r2] = regres3(t1rvll(idv,l1),t1rvll(idv,l2),1);
      y = b(1)*mnmx+b(2);
      plot(mnmx,y,'-','Color',colr(kd,:),'LineWidth',1);
      stxt = ['Slope = ' sprintf('%6.4f',b(1))]; % Slope text 
      ctxt = ['R^2 = ' sprintf('%6.4f',r2)];     % R^2 text
      text(mnmx(kd),y(kd),{stxt; ctxt},'Color',colr(kd,:), ...
           'FontSize',11,'FontWeight','bold','HorizontalAlignment', ...
           'center');
%
   end
%
   legend(hl,lld,'Location','southeast');
   axlim = axis;
   axis([0 axlim(2) 0 axlim(4)]);
   axis equal;
   xlabel('Visit 1 (ms)','FontSize',12,'FontWeight','bold');
   ylabel('Visit 2 (ms)','FontSize',12,'FontWeight','bold');
   ttxt = {'T1\rho '; [lleg{kl} ' Leg ']};
   title(ttxt,'FontSize',16,'FontWeight','bold');
%
   if kl==1
     print -dpsc2 -r600 -fillpage ...
           Results\NACOB_Final4\T1r_visit_corr4.ps;
   else
     print -dpsc2 -r600 -fillpage -append ...
           Results\NACOB_Final4\T1r_visit_corr4.ps;
   end
%
end
%
% Plots by Leg, Load, Compartment, and Bone
%
for kl = 1:2            % Leg
%
   ll = 8*kl;
%
   for kd = 1:2         % Load
%
       ld = 4*kd;
       
      for kc = 1:2      % Compartment
%
         lc = 2*kc;
%
         figure;
         orient landscape;
         hold on;
%
         for kb = 1:2   % Bone
%
            l1 = ll+ld+lc+kb-14;
            l2 = l1+16;
%
            if imrk
              text(t1rvllcb(:,l1),t1rvllcb(:,l2),subjt, ...
                   'Color',colr(kb,:), 'VerticalAlignment','middle', ...
                   'HorizontalAlignment','center');
            else
              hl(kb) = plot(t1rvllcb(:,l1),t1rvllcb(:,l2),'.', ...
                            'Color',colr(kb,:),'MarkerSize',8, ...
                            'LineWidth',1);
            end
            idv = ~isnan(t1rvllcb(:,l1));   % NaNs only on Visit 1
            [b,r2] = regres3(t1rvllcb(idv,l1),t1rvllcb(idv,l2),1);
            y = b(1)*mnmx+b(2);
            plot(mnmx,y,'-','Color',colr(kb,:),'LineWidth',1);
            stxt = ['Slope = ' sprintf('%6.4f',b(1))];     % Slope text 
            ctxt = ['R^2 = ' sprintf('%6.4f',r2)];         % R^2 text
            text(mnmx(kb),y(kb),{stxt; ctxt},'Color',colr(kb,:), ...
                 'FontSize',11,'FontWeight','bold', ...
                 'HorizontalAlignment','center');
%
         end            % End of bone loop (kb)
%
         if imrk
           legend(lbone,'Location','southeast');
         else
           legend(hl,lbone,'Location','southeast');
         end
         axlim = axis;
         axis([0 axlim(2) 0 axlim(4)]);
         axis equal;
         xlabel('Visit 1 (ms)','FontSize',12,'FontWeight','bold');
         ylabel('Visit 2 (ms)','FontSize',12,'FontWeight','bold');
         ttxt = {['T1\rho ' lleg{kl} ' Leg ']; [lld{kd} ' ' ...
                 lcmprt{kc}]};
         title(ttxt,'FontSize',16,'FontWeight','bold');
%
         print -dpsc2 -r600 -fillpage -append ...
               Results\NACOB_Final4\T1r_visit_corr4.ps;
%
      end               % End of compartment loop (kc)
   end                  % End of load loop (kd)
end                     % End of leg loop (kl)
%
% T2* Means
%
% mri_fitr4stat.xlsx has results from poor registrations set to -1
%
num = xlsread(fullfile(rdir,'mri_fitr4stat.xlsx'),'T2star');
%
% Indices
%
ivis = logical(num(:,2));              % Visit
ileg = logical(num(:,4));              % Leg
ild = logical(num(:,5));               % Load
icmprt = logical(num(:,6));            % Compartment
ibone = logical(num(:,7));             % Bone
t2s = num(:,12);                       % T1rho
%
% Change Invalid T1rho from -1 to NaN
%
t2s(t2s==-1) = NaN;
%
% Get Indices for Visits 1 and 2
%
ilegv(:,2) = ileg(ivis);
ilegv(:,1) = ileg(~ivis);
ildv(:,2) = ild(ivis);
ildv(:,1) = ild(~ivis);
icmprtv(:,2) = icmprt(ivis);
icmprtv(:,1) = icmprt(~ivis);
ibonev(:,2) = ibone(ivis);
ibonev(:,1) = ibone(~ivis);
t2sv(:,2) = t2s(ivis);
t2sv(:,1) = t2s(~ivis);
%
% Get Indices for Left and Right Legs
%
% Indices:
% 1 - Visit 1 Left Leg
% 2 - Visit 1 Right Leg
% 3 - Visit 2 Left Leg
% 4 - Visit 2 Right Leg
%
ildvl(:,4) = ildv(ilegv(:,2),2);
ildvl(:,3) = ildv(~ilegv(:,2),2);
ildvl(:,2) = ildv(ilegv(:,1),1);
ildvl(:,1) = ildv(~ilegv(:,1),1);
%
icmprtvl(:,4) = icmprtv(ilegv(:,2),2);
icmprtvl(:,3) = icmprtv(~ilegv(:,2),2);
icmprtvl(:,2) = icmprtv(ilegv(:,1),1);
icmprtvl(:,1) = icmprtv(~ilegv(:,1),1);
%
ibonevl(:,4) = ibonev(ilegv(:,2),2);
ibonevl(:,3) = ibonev(~ilegv(:,2),2);
ibonevl(:,2) = ibonev(ilegv(:,1),1);
ibonevl(:,1) = ibonev(~ilegv(:,1),1);
%
t2svl(:,4) = t2sv(ilegv(:,2),2);
t2svl(:,3) = t2sv(~ilegv(:,2),2);
t2svl(:,2) = t2sv(ilegv(:,1),1);
t2svl(:,1) = t2sv(~ilegv(:,1),1);
%
% Get Indices for Unloaded and Loaded
%
% Indices:
% 1 - Visit 1 Left Leg Unloaded
% 2 - Visit 1 Left Leg Loaded
% 3 - Visit 1 Right Leg Unloaded
% 4 - Visit 1 Right Leg Loaded
% 5 - Visit 2 Left Leg Unloaded
% 6 - Visit 2 Left Leg Loaded
% 7 - Visit 2 Right Leg Unloaded
% 8 - Visit 2 Right Leg Loaded
%
for k = 8:-2:2
   l = k/2;
   icmprtvll(:,k) = icmprtvl(ildvl(:,l),l);
   ibonevll(:,k) = ibonevl(ildvl(:,l),l);
   t2svll(:,k) = t2svl(ildvl(:,l),l);
end
%
for k = 7:-2:1
   l = (k+1)/2;
   icmprtvll(:,k) = icmprtvl(~ildvl(:,l),l);
   ibonevll(:,k) = ibonevl(~ildvl(:,l),l);
   t2svll(:,k) = t2svl(~ildvl(:,l),l);
end
%
% Get Indices for Lateral and Medial Compartments
%
% Indices:
% 1 - Visit 1 Left Leg Unloaded Lateral
% 2 - Visit 1 Left Leg Unloaded Medial
% 3 - Visit 1 Left Leg Loaded Lateral
% 4 - Visit 1 Left Leg Loaded Medial
% 5 - Visit 1 Right Leg Unloaded Lateral
% 6 - Visit 1 Right Leg Unloaded Medial
% 7 - Visit 1 Right Leg Loaded Lateral
% 8 - Visit 1 Right Leg Loaded Medial
% 9 - Visit 2 Left Leg Unloaded Lateral
% 10 - Visit 2 Left Leg Unloaded Medial
% 11 - Visit 2 Left Leg Loaded Lateral
% 12 - Visit 2 Left Leg Loaded Medial
% 13 - Visit 2 Right Leg Unloaded Lateral
% 14 - Visit 2 Right Leg Unloaded Medial
% 15 - Visit 2 Right Leg Loaded Lateral
% 16 - Visit 2 Right Leg Loaded Medial
%
for k = 16:-2:2
   l = k/2;
   ibonevllc(:,k) = ibonevll(icmprtvll(:,l),l);
   t2svllc(:,k) = t2svll(icmprtvll(:,l),l);
end
%
for k = 15:-2:1
   l = (k+1)/2;
   ibonevllc(:,k) = ibonevll(~icmprtvll(:,l),l);
   t2svllc(:,k) = t2svll(~icmprtvll(:,l),l);
end
%
% Get Indices for the Femur and Tibia
%
% Indices:
% 1 - Visit 1 Left Leg Unloaded Lateral Femur
% 2 - Visit 1 Left Leg Unloaded Lateral Tibia
% 3 - Visit 1 Left Leg Unloaded Medial Femur
% 4 - Visit 1 Left Leg Unloaded Medial Tibia
% 5 - Visit 1 Left Leg Loaded Lateral Femur
% 6 - Visit 1 Left Leg Loaded Lateral Tibia
% 7 - Visit 1 Left Leg Loaded Medial Femur
% 8 - Visit 1 Left Leg Loaded Medial Tibia
% 9 - Visit 1 Right Leg Unloaded Lateral Femur
% 10 - Visit 1 Right Leg Unloaded Lateral Tibia
% 11 - Visit 1 Right Leg Unloaded Medial Femur
% 12 - Visit 1 Right Leg Unloaded Medial Tibia
% 13 - Visit 1 Right Leg Loaded Lateral Femur
% 14 - Visit 1 Right Leg Loaded Lateral Tibia
% 15 - Visit 1 Right Leg Loaded Medial Femur
% 16 - Visit 1 Right Leg Loaded Medial Tibia
% 17 - Visit 2 Left Leg Unloaded Lateral Femur
% 18 - Visit 2 Left Leg Unloaded Lateral Tibia
% 19 - Visit 2 Left Leg Unloaded Medial Femur
% 20 - Visit 2 Left Leg Unloaded Medial Tibia
% 21 - Visit 2 Left Leg Loaded Lateral Femur
% 22 - Visit 2 Left Leg Loaded Lateral Tibia
% 23 - Visit 2 Left Leg Loaded Medial Femur
% 24 - Visit 2 Left Leg Loaded Medial Tibia
% 25 - Visit 2 Right Leg Unloaded Lateral Femur
% 26 - Visit 2 Right Leg Unloaded Lateral Tibia
% 27 - Visit 2 Right Leg Unloaded Medial Femur
% 28 - Visit 2 Right Leg Unloaded Medial Tibia
% 29 - Visit 2 Right Leg Loaded Lateral Femur
% 30 - Visit 2 Right Leg Loaded Lateral Tibia
% 31 - Visit 2 Right Leg Loaded Medial Femur
% 32 - Visit 2 Right Leg Loaded Medial Tibia
%
for k = 32:-2:2
   l = k/2;
   t2svllcb(:,k) = t2svllc(ibonevllc(:,l),l);
end
%
for k = 31:-2:1
   l = (k+1)/2;
   t2svllcb(:,k) = t2svllc(~ibonevllc(:,l),l);
end
%
% Plot Leg and Load
%
mnmx = [min(min(t2svll)) max(max(t2svll))];
%
for kl = 1:2            % Leg
%
   l = 2*kl;
   figure;
   orient landscape;
   hold on;
%
   for kd = 1:2         % Load
%
      l1 = l+kd-2;
      l2 = l1+4;
%
      hl(kd) = plot(t2svll(:,l1),t2svll(:,l2),'.','Color', ...
                    colr(kd,:),'MarkerSize',8,'LineWidth',1);
      idv = ~isnan(t2svll(:,l2));      % NaNs only on Visit 2
      [b,r2] = regres3(t2svll(idv,l1),t2svll(idv,l2),1);
      y = b(1)*mnmx+b(2);
      plot(mnmx,y,'-','Color',colr(kd,:),'LineWidth',1);
      stxt = ['Slope = ' sprintf('%6.4f',b(1))]; % Slope text 
      ctxt = ['R^2 = ' sprintf('%6.4f',r2)];     % R^2 text
      text(mnmx(kd),y(kd),{stxt; ctxt},'Color',colr(kd,:), ...
           'FontSize',11,'FontWeight','bold','HorizontalAlignment', ...
           'center');
%
   end
%
   legend(hl,lld,'Location','southeast');
   axlim = axis;
   axis([0 axlim(2) 0 axlim(4)]);
   axis equal;
   xlabel('Visit 1 (ms)','FontSize',12,'FontWeight','bold');
   ylabel('Visit 2 (ms)','FontSize',12,'FontWeight','bold');
   ttxt = {'T2* '; [lleg{kl} ' Leg ']};
   title(ttxt,'FontSize',16,'FontWeight','bold');
%
   if kl==1
     print -dpsc2 -r600 -fillpage ...
           Results\NACOB_Final4\T2s_visit_corr4.ps;
   else
     print -dpsc2 -r600 -fillpage -append ...
           Results\NACOB_Final4\T2s_visit_corr4.ps;
   end
%
end
%
% Plots by Leg, Load, Compartment, and Bone
%
for kl = 1:2            % Leg
%
   ll = 8*kl;
%
   for kd = 1:2         % Load
%
       ld = 4*kd;
       
      for kc = 1:2      % Compartment
%
         lc = 2*kc;
%
         figure;
         orient landscape;
         hold on;
%
         for kb = 1:2   % Bone
%
            l1 = ll+ld+lc+kb-14;
            l2 = l1+16;
%
            if imrk
              text(t2svllcb(:,l1),t2svllcb(:,l2),subjt, ...
                   'Color',colr(kb,:), 'VerticalAlignment','middle', ...
                   'HorizontalAlignment','center');
            else
              hl(kb) = plot(t2svllcb(:,l1),t2svllcb(:,l2),'.', ...
                            'Color',colr(kb,:),'MarkerSize',8, ...
                            'LineWidth',1);
            end
            idv = ~isnan(t2svllcb(:,l2));   % NaNs only on Visit 2
            [b,r2] = regres3(t2svllcb(idv,l1),t2svllcb(idv,l2),1);
            y = b(1)*mnmx+b(2);
            plot(mnmx,y,'-','Color',colr(kb,:),'LineWidth',1);
            stxt = ['Slope = ' sprintf('%6.4f',b(1))];     % Slope text 
            ctxt = ['R^2 = ' sprintf('%6.4f',r2)];         % R^2 text
            text(mnmx(kb),y(kb),{stxt; ctxt},'Color',colr(kb,:), ...
                 'FontSize',11,'FontWeight','bold', ...
                 'HorizontalAlignment','center');
%
         end            % End of bone loop (kb)
%
         if imrk
           legend(lbone,'Location','southeast');
         else
           legend(hl,lbone,'Location','southeast');
         end
         axlim = axis;
         axis([0 axlim(2) 0 axlim(4)]);
         axis equal;
         xlabel('Visit 1 (ms)','FontSize',12,'FontWeight','bold');
         ylabel('Visit 2 (ms)','FontSize',12,'FontWeight','bold');
         ttxt = {['T2* ' lleg{kl} ' Leg ']; [lld{kd} ' ' lcmprt{kc}]};
         title(ttxt,'FontSize',16,'FontWeight','bold');
%
         print -dpsc2 -r600 -fillpage -append ...
               Results\NACOB_Final4\T2s_visit_corr4.ps;
%
      end               % End of compartment loop (kc)
   end                  % End of load loop (kd)
end                     % End of leg loop (kl)
%
return