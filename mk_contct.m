function [tptsc,fptsc] = mk_contct(tpts,fpts,ipts,idcs,iplt,plt_title)
%MK_CONTCT Finds the distances between a pair of piece-wise linear
%          lines and moves the lines to be in contact midway between
%          the lines.
%
%          [TPTSC,FPTSC] = MK_CONTCT(TPTS,FPTS,IPTS,IDCS) Given a matrix
%          with the two- or three-dimensional coordinates (in columns)
%          of the points forming a piece-wise linear line, TPTS, a
%          similar matrix forming a second piece-wise linear line, FPTS,
%          coordinates (in columns) of the intersection points between
%          the two lines, IPTS, and a two column matrix with the index
%          to the segments within the piece-wise linear line with an
%          intersection (first line in column 1 and second line in
%          column 2), returns the two piece-wise linear lines with
%          additional points and regions of overlap in contact, TPTS,
%          and FPTS.
%
%          [TPTSC,FPTSC] = MK_CONTCT(TPTS,FPTS,IPTS,IDCS,IPLT) If IPLT
%          is true (or nonzero), plots the original lines, intersection
%          points, closest points between the lines, and the updated
%          lines in contact.
%
%          [TPTSC,FPTSC] = MK_CONTCT(TPTS,FPTS,IPTS,IDCS,IPLT,PLT_TITLE)
%          Given a plot tile, PLT_TITLE, adds a title to the plot.
%
%          NOTES:  1.  The points along the lines are sorted based on
%                  the X coordinates in 2-D and Y coordinates in 3-D.
%
%                  2.  The M-files dis2lin.m and pts2lin.m must be in
%                  the current path or directory.
%
%                  3.  The function works with either two- or three-
%                  dimensional coordinates.
%
%          05-Apr-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  error(' *** ERROR in MK_CONTCT:  Four inputs are required!');
end
%
if (nargin<5)
  iplt = false;
end
%
if (nargin<6)
  plt_title = [];
end
%
% Initialize Arrays
%
dom = size(tpts,2);     % Size of domain
if dom==2
  isrt = 1;             % Sort 2-D points along the line by X (1st column)
  dom2 = true;
else
  isrt = 2;             % Sort 3-D points along the line by Y (2nd column)
  dom2 = false;
end
%
% Plot Lines
%
if iplt
%
    figure;
    orient landscape;
%
  if dom2
%
    plot(fpts(:,1),fpts(:,2),'b.-','LineWidth',1.5,'MarkerSize',12);
    hold on;
    plot(tpts(:,1),tpts(:,2),'g.-','LineWidth',1.5,'MarkerSize',12);
    plot(ipts(:,1),ipts(:,2),'ro','LineWidth',1.5,'MarkerSize',8);
    set(gca,'XDir','reverse');
    set(gca,'YDir','reverse');
    view(2);
%
  else
%
    plot3(fpts(:,1),fpts(:,2),fpts(:,3),'b.-','LineWidth', ...
          1.5,'MarkerSize',12);
    hold on;
    plot3(tpts(:,1),tpts(:,2),tpts(:,3),'g.-','LineWidth', ...
          1.5,'MarkerSize',12);
    plot3(ipts(:,1),ipts(:,2),ipts(:,3),'ro','LineWidth',1.5, ...
          'MarkerSize',8);
    view(-90,0);
%
  end
%
  axis equal;
  title(plt_title,'FontSize',16,'FontWeight','bold');
%
end
%
% Loop through Pairs of Intersection Points
%
fptsc = fpts;
tptsc = tpts;
%
ni = size(ipts,1);      % Number of intersecting points
%
for km = 1:ni/2
   mi = 2*km;
   idfx = fpts(:,isrt)>ipts(mi-1,isrt)&fpts(:,isrt)<ipts(mi,isrt);
   idtx = tpts(:,isrt)>ipts(mi-1,isrt)&tpts(:,isrt)<ipts(mi,isrt);
   idfxc = fptsc(:,isrt)>ipts(mi-1,isrt)&fptsc(:,isrt)<ipts(mi,isrt);
   idtxc = tptsc(:,isrt)>ipts(mi-1,isrt)&tptsc(:,isrt)<ipts(mi,isrt);
   nf = nnz(idfx);      % Number of points on the femur line in the overlap
   nt = nnz(idtx);      % Number of points on the tibia line in the overlap
%
   idcs(mi-1:mi,:) = sort(idcs(mi-1:mi,:)); % Make sure both lines are going in the same direction
%
% Overlapping Femur Points
%
   if nf>0
%
% Find Shortest Distance Between the Lines
%
     idtl = idcs(mi-1,1):idcs(mi,1)+1;
     txyz = tpts(idtl,:);
     fxyz = fpts(idfx,:);
     [~,idf,xyzf] = dis2lin(txyz,fxyz);
%
% Find Midpoints and Update Coordinates
%
     fmdpts = (fxyz(idf,:)+xyzf)/2;    % Midpoints
     idfxc(idfxc) = idf;
     fptsc(idfxc,:) = fmdpts;          % Update femur points positions
%
% Plot Femur Midpoints
%
     if iplt
       if dom2
         plot([fxyz(idf,1)'; xyzf(:,1)'],[fxyz(idf,2)'; xyzf(:,2)'], ...
              'm-','LineWidth',1.5);
         plot(fmdpts(:,1)',fmdpts(:,2)','rs','LineWidth',1.5, ...
              'MarkerSize',8);
       else
         plot3([fxyz(idf,1)'; xyzf(:,1)'],[fxyz(idf,2)'; ...
               xyzf(:,2)'],[fxyz(idf,3)'; xyzf(:,3)'],'m-', ...
               'LineWidth',1.5);
         plot3(fmdpts(:,1)',fmdpts(:,2)',fmdpts(:,3)','rs', ...
               'LineWidth',1.5,'MarkerSize',8);
       end
     end
%
   end
%
% Overlapping Tibia Points
%
   if nt>0
%
% Find Shortest Distance Between the Lines
%
     idfl = idcs(mi-1,2):idcs(mi,2)+1;
     txyz = tpts(idtx,:);
     fxyz = fpts(idfl,:);
     [~,idt,xyzt] = dis2lin(fxyz,txyz);
%
% Find Midpoints and Update Coordinates
%
     tmdpts = (txyz(idt,:)+xyzt)/2;    % Midpoints
     idtxc(idtxc) = idt;
     tptsc(idtxc,:) = tmdpts;          % Update tibia points positions
%
% Plot Tibia Midpoints
%
     if iplt
       if dom2
         plot([txyz(idt,1)'; xyzt(:,1)'],[txyz(idt,2)'; xyzt(:,2)'], ...
              'c-','LineWidth',1.5);
         plot(tmdpts(:,1),tmdpts(:,2),'rs','LineWidth',1.5, ...
              'MarkerSize',8);
       else
         plot3([txyz(idt,1)'; xyzt(:,1)'],[txyz(idt,2)'; ...
               xyzt(:,2)'],[txyz(idt,3)'; xyzt(:,3)'],'c-', ...
               'LineWidth',1.5);
         plot3(tmdpts(:,1),tmdpts(:,2),tmdpts(:,3),'rs', ...
               'LineWidth',1.5,'MarkerSize',8);
       end
     end
%
   end
%
% Add in New Points
%
   if nf>0
     if all(~idtxc)     % No tibia points in overlap
       mdpt = (ipts(mi-1,:)+ipts(mi,:))/2;
       ntp = size(tptsc,1);
       d = tptsc-repmat(mdpt,ntp,1);
       d = sum(d.*d,2);
       [~,idtr] = min(d);
     else
       idtr = find(idtxc);
     end
     idtr = min(idtr)-1:max(idtr)+1;
     txyz = [tptsc(idtr,:); fmdpts];
     txyz = sortrows(txyz,sign(txyz(2,isrt)-txyz(1,isrt))*isrt);
     tptsc = [tptsc(1:idtr(1)-1,:); txyz; tptsc(idtr(end)+1:end,:)];
   end
%
   if nt>0
     if all(~idfxc)     % No femur points in overlap
       mdpt = (ipts(mi-1,:)+ipts(mi,:))/2;
       nfp = size(fptsc,1);
       d = fptsc-repmat(mdpt,nfp,1);
       d = sum(d.*d,2);
       [~,idfr] = min(d);
     else
       idfr = find(idfxc);
     end
     idfr = min(idfr)-1:max(idfr)+1;
     fxyz = [fptsc(idfr,:); tmdpts];
     fxyz = sortrows(fxyz,sign(fxyz(2,isrt)-fxyz(1,isrt))*isrt);
     fptsc = [fptsc(1:idfr(1)-1,:); fxyz; fptsc(idfr(end)+1:end,:)];
   end
%
end
%
% Plot New Lines in Contact
%
if iplt
%
  if dom2
%
    plot(fptsc(:,1),fptsc(:,2),'bv:','LineWidth',1.5,'MarkerSize',8);
    plot(tptsc(:,1),tptsc(:,2),'g^:','LineWidth',1.5,'MarkerSize',8);
%
  else
%
    plot3(fptsc(:,1),fptsc(:,2),fptsc(:,3),'bv:','LineWidth', ...
          1.5,'MarkerSize',8);
    plot3(tptsc(:,1),tptsc(:,2),tptsc(:,3),'g^:','LineWidth', ...
          1.5,'MarkerSize',8);
    axis tight;
  end
%
end
%
return