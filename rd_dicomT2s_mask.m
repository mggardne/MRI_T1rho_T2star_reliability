%#######################################################################
%
%                * ReaD DICOM T2* Volume Data Program *
%
%          M-File which reads the DICOM images from a T2* series.  The
%     image data is put into volume matrices and the volumes from
%     different echo times are registered to the echo time = 5 ms 
%     volume for T2*.  The registered volume data and series
%     information are saved to a MAT file.
%
%          A single T2* series of echo times must be manually selected.
%
%     NOTES:  1.  For use with subject 08-AS on Visit 2 for the T2* for
%             the unloaded left leg.
%
%             2.  Matlab MAT file dicom_lst2.mat must be in the current
%             directory or path.
%
%             3.  elastix.exe must be installed and executable on this
%             computer.
%
%             4.  The Matlab path must include the following two paths:
% C:\Users\mggardne\BRUCE\Risk_Fac\CACL\MRI_data\MelastiX\yamlmatlab-master
% C:\Users\mggardne\BRUCE\Risk_Fac\CACL\MRI_data\MelastiX\matlab_elastix-master\code
%             These are required to run elastix using Rob Campbell's
%             MelastiX.
%
%             5.  Elastix parameter file, Parameters_RigidBody.txt,
%             must in the current directory.
%
%             6.  The M-files dftreg.m and dftregistration.m must be in
%             the current path or directory.
%
%             7.  The Matlab Image Processing Toolbox is required.
%
%             8.  The echo times must be certain values and in
%             ascending order.
%
%     03-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
rad2deg = 180/pi;       % Radians to degrees
%
% Load Image File Data from dicom_lst2.mat
%
load dicom_lst2.mat afiles ddirs etn ets isz nimages pspc sn stxt;
%
% Use Series Descriptions to Find T2* Series
%
ids = contains(stxt,'*')|contains(stxt,'T2','IgnoreCase',true); % T2*
n = sum(ids);           % Number of possible series times to analyze
%
sn = sn(ids);           % Series numbers
etss = strrep(ets(ids),'(1)',' ms');   % T2* echo times
stxts = stxt(ids);                     % T2* series
stxtss = cell(n,1);
for k = 1:n
   stxtss{k} = ['Series ' int2str(sn(k)) ', ' stxts{k}, ...
                ', echo time = ' etss{k}];
end
%
[idse,ok] = listdlg('ListString',stxtss,'SelectionMode','multiple', ...
                    'Name','T2* Series','PromptString', ...
                    'Select a T2* Series','ListSize',[300 400]);
%
if ok<1
  warning(' *** WARNING in rd_dicomT2*:  No series selected!');
  return
end
%
idse = idse';
ids = find(ids);
ids = ids(idse);        % Index to series to analyze
netn = size(ids,1);     % Number of echo times to analyze
%
etns = [etn{ids}]';     % T2* echo times
[etns,idss] = sort(etns);              % Sorted echo times
ids = ids(idss);        % Sorted index to T2* series
idse = idse(idss);      % Sorted index within T2* series
%
ns = 1;                 % Number of sets of T2* sequences
%
% Get Series Index and Echo Times as Matrices
%
etnm = cell2mat(etn(ids));             % Echo times for series
%
[~,idsrt] = sort(etnm);
if all(idsrt~=(1:netn)')
  error(' *** rd_dicom:  Echo times are not in order!');
end
%
% Get Index to Echo Time = 5 ms
%
d = etnm-5;
[~,id5] = min(d.*d);
id5 = id5';
%
% Get Mask Parameters
%
idmx = 54:288;
idmy = 49:255;
% idmx = 60:244;
% idmy = 65:278;
iszsp = [size(idmx,2) size(idmy,2)];
%
% Get Additional T2* Variables
%
nfiles = nimages(ids(:,1));  % Numbers of T2* files
iszss = isz(ids(:,1),:);     % Image sizes in pixels
pspcss = pspc(ids(:,1),:);   % Pixel sizes
%
ddirss = ddirs(ids);         % Subdirectories for T2* series
afiless = afiles(ids);       % T2* files
snss = sn(idse);             % Series numbers
stxts = stxts(idse);         % T2* series
%
% Loop through the T2* Series
%
nreg = ns*(netn-1);     % Number of volume registrations
tms = zeros(nreg,10);   % 3D and 2D registration translations/rotations
ett = cell(nreg,1);     % Registration echo times
ss = cell(nreg,1);      % Series
fits = string(repmat('T2star',nreg,1));% Type of fit (T1rho or T2star)
%
for k = 1:ns
%
% Get Series Specific Indices
%
   ddir = ddirss(:,k);  % Subdirectories for series
%
   nfile = nfiles(k);   % Number of image files (slices) in this series
   etns = etnm(:,k);    % T2* echo times
   iszs = iszss(k,:);   % Image size
   pspcs = pspcss(k,:); % Pixel size
%
   sns = snss(:,k);     % Series numbers
   snt = int2str(sns(1));              % First series number as text
   psfile = ['S' snt '.ps'];           % PS file name
%
% Get Images and Maximum Scaled Image Values
%
   valmx = zeros(nfile*netn,1);
   v = zeros([iszs nfile netn]);
%
   for l = 1:netn
%
% Get Image File Names
%
      fnams = afiless{l,k};            % T2* files
%
      for m = 1:nfile
         info = dicominfo(fullfile(ddir{l},fnams{m}));
         sl = double(info.RescaleSlope);
         y0 = double(info.RescaleIntercept);
         img = dicomread(info);
         img = sl*double(img)+y0;
         v(:,:,m,l) = img;
         idx = l*nfile+m-nfile;
         valmx(idx) = max(img(:));     % Slice maximum
      end
   end
%
%    vt = v(101:200,91:250,:,:);
%    iszsp = [100 160];
%    vt = v(49:255,54:288,:,:);
%    iszsp = [207 235];
   vt = v(idmx,idmy,:,:);
%    iszsp = [235 207];
%
   scmx = 10*fix(max(valmx)/10);       % Round maximum value down
%
% Register the Different Echo Time Images to the 5 ms Echo Time Image
%
   rsl(2,1) = floor(nfile/4);
   rsl = [rsl(2,1); 3*rsl(2,1)];       % Slices for 2D registration
%
   pxmx = 3*iszsp(2);                  % Maximum X pixels
   px = 0:23:pxmx;                     % Grid of 23 pixels
   nx = size(px,2);                    % Number of X grid lines
   pymx = iszsp(1);                    % Maximum Y pixels
   py = 0:47:pymx;                     % Grid of 47 pixels
   ny = size(py,2);                    % Number of Y grid lines
%
   t2 = cell(2,1);      % 2D translations
%
% Register Echo Times Before 5 ms Echo Times
%
   idsrb = id5(k):-1:1;      % Echo times before 5 ms echo time
   nb = size(idsrb,2)-1;     % Number of before echo times to register
%
   for l = 1:nb
%
      n1 = idsrb(l);    % Index to fixed image
      n2 = idsrb(l+1);  % Index to moving image
%
      o = netn*k-netn-k+l+1;
      ss{o} = ['S' snt];
      ett{o} = ['Echo time ' sprintf('%.2f',etns(n2)), ...
                ' to ' sprintf('%.2f',etns(n1)) ' ms'];
%
      [regt,t] = elastix(vt(:,:,:,n2),vt(:,:,:,n1),[], ...
                         'Parameters_RigidBody.txt');      % 3D
%
      s = t.TransformParameters{1};
      s.Size(1:2) = iszs;
      s.CenterOfRotationPoint(1:2) = s.CenterOfRotationPoint(1:2)+ ...
                                     [idmx(1) idmy(1)];
      t.TransformParameters{1} = s;
%
      reg = transformix(v(:,:,:,n2),t);     % Transform full image matrix
%
      t = t.TransformParameters{1};

      rx = t.TransformParameters(1);
      ry = t.TransformParameters(2);
      rz = t.TransformParameters(3);
      tx = t.TransformParameters(4);
      ty = t.TransformParameters(5);
      tz = t.TransformParameters(6);
      tn = [tx; ty; tz; rad2deg*[rx; ry; rz]];
      t = sprintf(['tx = %.1f, ty = %.1f, tz = %.1f, rx =  %.1f, ', ...
                   'ry =  %.1f, rz =  %.1f'],tn);
%
      [tf1,dr1,dc1] = dftreg(vt(:,:,rsl(1),n1), ...
                             vt(:,:,rsl(1),n2),100);       % 2D
      t2{1} = sprintf('tx = %.1f, ty = %.1f',dc1,dr1);
      [tf2,dr2,dc2] = dftreg(vt(:,:,rsl(2),n1), ...
                             vt(:,:,rsl(2),n2),100);       % 2D
      t2{2} = sprintf('tx = %.1f, ty = %.1f',dc2,dr2);
%
      tms(o,:) = [tn',dc1,dr1,dc2,dr2];
%
% Plot Registration
%
      for m = 1:2
%
% Plots with Grids
%
        vp = squeeze(vt(:,:,rsl(m),n1:-1:n2));
        irng = [min(vp(:)) max(vp(:))];
        figure;
        subplot(2,1,1);
        montage({vp(:,:,1),vp(:,:,2),regt(:,:,rsl(m))}, ...
                'Size',[1 3],'DisplayRange',irng,'ThumbnailSize',iszsp);
        orient landscape;
        hold on;
        plot3(repmat(px,2,1),repmat([0; pymx],1,nx),ones(2,nx),'r-');
        plot3(repmat([0; pxmx],1,ny),repmat(py,2,1),ones(2,ny),'r-');
        xlabel(t,'FontSize',12,'FontWeight','bold');
        ylabel('3D','FontSize',12,'FontWeight','bold');
        title({['Series ' snt]; ett{o}; ['Slice ' int2str(rsl(m))]}, ...
              'FontSize',16,'FontWeight','bold');
%
        subplot(2,1,2);
        if m==1
          montage({vp(:,:,1),vp(:,:,2),tf1},'Size', ...
                  [1 3],'DisplayRange',irng,'ThumbnailSize',iszsp);
        else
          montage({vp(:,:,1),vp(:,:,2),tf2},'Size', ...
                  [1 3],'DisplayRange',irng,'ThumbnailSize',iszsp);
        end
        colormap gray;
        n3 = nb-l+1;
        brighten(n3*n3/32);
        hold on;
        plot3(repmat(px,2,1),repmat([0; pymx],1,nx),ones(2,nx),'r-');
        plot3(repmat([0; pxmx],1,ny),repmat(py,2,1),ones(2,ny),'r-');
        xlabel(t2{m},'FontSize',12,'FontWeight','bold');
        ylabel('2D','FontSize',12,'FontWeight','bold');
%
        if l==1&&m==1
          print('-dpsc2','-r600','-fillpage',psfile);
        else
          print('-dpsc2','-r600','-fillpage','-append',psfile);
        end
%
% Plots Showing Differences as Color
%
        if m==1
          hf1 = figure;
          orient landscape;
          subplot(2,3,1);
          imshowpair(vp(:,:,1),vp(:,:,2));
          title({'Unregistered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,2);
          imshowpair(vp(:,:,1),regt(:,:,rsl(m)));
          xlabel(t,'FontSize',11,'FontWeight','bold');
          title({'3D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,3);
          imshowpair(vp(:,:,1),tf1);
          xlabel(t2{m},'FontSize',11,'FontWeight','bold');
          title({'2D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
        else
          figure(hf1);
          subplot(2,3,4);
          imshowpair(vp(:,:,1),vp(:,:,2));
          title({'Unregistered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,5);
          imshowpair(vp(:,:,1),regt(:,:,rsl(m)));
          xlabel(t,'FontSize',11,'FontWeight','bold');
          title({'3D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,6);
          imshowpair(vp(:,:,1),tf2);
          xlabel(t2{m},'FontSize',11,'FontWeight','bold');
          title({'2D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          sgtitle({['Series ' snt]; ett{o}}, ...
                  'FontSize',16,'FontWeight','bold');
%
          print('-dpsc2','-r600','-fillpage','-append',psfile);
%
        end
%
% Plots Full Image Compared to Masked Image Showing Differences as Color
%
        figure;
        orient landscape;
        subplot(1,3,1);
        imshowpair(v(:,:,rsl(m),n1),reg(:,:,rsl(m)));
        title({'Full Images'; ['Slice ', ...
              int2str(rsl(m))]},'FontSize',12,'FontWeight','bold');
%
        subplot(1,3,2);
        imshowpair(vt(:,:,rsl(m),n1),regt(:,:,rsl(m)));
        title({'Masked Images'; ['Slice ', ...
              int2str(rsl(m))]},'FontSize',12,'FontWeight','bold');
%
        subplot(1,3,3);
        imshowpair(reg(idmx,idmy,rsl(m)),regt(:,:,rsl(m)));
        title({'Full Compared to Masked'; ['Slice ', ...
              int2str(rsl(m))]},'FontSize',12,'FontWeight','bold');
%
        sgtitle({['Series ' snt]; ett{o}}, ...
                'FontSize',16,'FontWeight','bold');
%
        print('-dpsc2','-r600','-fillpage','-append',psfile);
%
     end
%
     v(:,:,:,n2) = reg;  % Replace moving image with registered image
%
   end
%
% Register Echo Times After 5 ms Echo Times
%
   idsra = id5(k):1:netn;    % Echo times after 5 ms echo time
   na = size(idsra,2)-1;     % Number of after echo times to register
%
   for l = 1:na
%
      n1 = idsra(l);    % Index to fixed image
      n2 = idsra(l+1);  % Index to moving image
%
      o = netn*k-netn-k+l+nb+1;
      ss{o} = ['S' snt];
      ett{o} = ['Echo time ' sprintf('%.2f',etns(n2)), ...
                ' to ' sprintf('%.2f',etns(n1)) ' ms'];
%
      [regt,t] = elastix(vt(:,:,:,n2),vt(:,:,:,n1),[], ...
                         'Parameters_RigidBody.txt');      % 3D
%
      s = t.TransformParameters{1};
      s.Size(1:2) = iszs;
      s.CenterOfRotationPoint(1:2) = s.CenterOfRotationPoint(1:2)+ ...
                                     [idmx(1) idmy(1)];
      t.TransformParameters{1} = s;
%
      reg = transformix(v(:,:,:,n2),t);     % Transform full image matrix
%
      t = t.TransformParameters{1};

      rx = t.TransformParameters(1);
      ry = t.TransformParameters(2);
      rz = t.TransformParameters(3);
      tx = t.TransformParameters(4);
      ty = t.TransformParameters(5);
      tz = t.TransformParameters(6);
      tn = [tx; ty; tz; rad2deg*[rx; ry; rz]];
      t = sprintf(['tx = %.1f, ty = %.1f, tz = %.1f, rx =  %.1f, ', ...
                   'ry =  %.1f, rz =  %.1f'],tn);
%
      [tf1,dr1,dc1] = dftreg(vt(:,:,rsl(1),n1), ...
                             vt(:,:,rsl(1),n2),100);       % 2D
      t2{1} = sprintf('tx = %.1f, ty = %.1f',dc1,dr1);
      [tf2,dr2,dc2] = dftreg(vt(:,:,rsl(2),n1), ...
                             vt(:,:,rsl(2),n2),100);       % 2D
      t2{2} = sprintf('tx = %.1f, ty = %.1f',dc2,dr2);
%
      tms(o,:) = [tn',dc1,dr1,dc2,dr2];
%
% Plot Registration
%
      for m = 1:2
%
% Plots with Grids
%
        vp = squeeze(vt(:,:,rsl(m),n1:n2));
        irng = [min(vp(:)) max(vp(:))];
        figure;
        subplot(2,1,1);
        montage({vp(:,:,1),vp(:,:,2),regt(:,:,rsl(m))}, ...
                'Size',[1 3],'DisplayRange',irng,'ThumbnailSize',iszsp);
        orient landscape;
        hold on;
        plot3(repmat(px,2,1),repmat([0; pymx],1,nx),ones(2,nx),'r-');
        plot3(repmat([0; pxmx],1,ny),repmat(py,2,1),ones(2,ny),'r-');
        xlabel(t,'FontSize',12,'FontWeight','bold');
        ylabel('3D','FontSize',12,'FontWeight','bold');
        title({['Series ' snt]; ett{o}; ['Slice ' int2str(rsl(m))]}, ...
              'FontSize',16,'FontWeight','bold');
%
        subplot(2,1,2);
        if m==1
          montage({vp(:,:,1),vp(:,:,2),tf1},'Size', ...
                  [1 3],'DisplayRange',irng,'ThumbnailSize',iszsp);
        else
          montage({vp(:,:,1),vp(:,:,2),tf2},'Size', ...
                  [1 3],'DisplayRange',irng,'ThumbnailSize',iszsp);
        end
        colormap gray;
        n3 = l+2;
        brighten(n3*n3/32);
        hold on;
        plot3(repmat(px,2,1),repmat([0; pymx],1,nx),ones(2,nx),'r-');
        plot3(repmat([0; pxmx],1,ny),repmat(py,2,1),ones(2,ny),'r-');
        xlabel(t2{m},'FontSize',12,'FontWeight','bold');
        ylabel('2D','FontSize',12,'FontWeight','bold');
%
        print('-dpsc2','-r600','-fillpage','-append',psfile);
%
% Plots Showing Differences as Color
%
        if m==1
          hf1 = figure;
          orient landscape;
          subplot(2,3,1);
          imshowpair(vp(:,:,1),vp(:,:,2));
          title({'Unregistered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,2);
          imshowpair(vp(:,:,1),regt(:,:,rsl(m)));
          xlabel(t,'FontSize',11,'FontWeight','bold');
          title({'3D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,3);
          imshowpair(vp(:,:,1),tf1);
          xlabel(t2{m},'FontSize',11,'FontWeight','bold');
          title({'2D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
        else
          figure(hf1);
          subplot(2,3,4);
          imshowpair(vp(:,:,1),vp(:,:,2));
          title({'Unregistered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,5);
          imshowpair(vp(:,:,1),regt(:,:,rsl(m)));
          xlabel(t,'FontSize',11,'FontWeight','bold');
          title({'3D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,6);
          imshowpair(vp(:,:,1),tf2);
          xlabel(t2{m},'FontSize',11,'FontWeight','bold');
          title({'2D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          sgtitle({['Series ' snt]; ett{o}}, ...
                  'FontSize',16,'FontWeight','bold');
%
          print('-dpsc2','-r600','-fillpage','-append',psfile);
%
        end
%
% Plots Full Image Compared to Masked Image Showing Differences as Color
%
        figure;
        orient landscape;
        subplot(1,3,1);
        imshowpair(v(:,:,rsl(m),n1),reg(:,:,rsl(m)));
        title({'Full Images'; ['Slice ', ...
              int2str(rsl(m))]},'FontSize',12,'FontWeight','bold');
%
        subplot(1,3,2);
        imshowpair(vt(:,:,rsl(m),n1),regt(:,:,rsl(m)));
        title({'Masked Images'; ['Slice ', ...
              int2str(rsl(m))]},'FontSize',12,'FontWeight','bold');
%
        subplot(1,3,3);
        imshowpair(reg(idmx,idmy,rsl(m)),regt(:,:,rsl(m)));
        title({'Full Compared to Masked'; ['Slice ', ...
              int2str(rsl(m))]},'FontSize',12,'FontWeight','bold');
%
        sgtitle({['Series ' snt]; ett{o}}, ...
                'FontSize',16,'FontWeight','bold');
%
        print('-dpsc2','-r600','-fillpage','-append',psfile);
%
      end
%
      v(:,:,:,n2) = reg;  % Replace moving image with registered image
%
   end
%
% Save Registered Image Data to a MAT File
%
   fnams = afiless{:,k};               % T2* files
   nsls = nfile;                       % Number of slices
   nfile = nfile*netn;                 % Total number of files
   st = stxts{1,k};                    % Series description
%
   matfile = ['T2star_S' snt '.mat'];
%
   save(matfile,'ddir','etns','fnams','id5','iszs','nfile','nsls', ...
        'netn','pspcs','scmx','sns','snt','st','v');
end
%
% Write Registration Transformations to the Spreadsheet
%
dirstr = split(pwd,filesep);
dirstr = [dirstr{end-1} '_' dirstr{end}];
dirstr = strrep(dirstr,' ','_');
%
xlsnam = [dirstr '.xlsx'];
%
colnam1 = {'Series','MRI_Type','Registered_MRI_Times'};
%
t1 = table(string(ss),fits,string(ett),'VariableNames',colnam1);
%
colnam2 = {'tx_3D','ty_3D','tz_3D','rx_3D','ry_3D','rz_3D', ...
           'tx_2Dslice1','ty_2Dslice1','tx_2Dslice2','ty_2Dslice2'};
%
t2 = array2table(tms,'VariableNames',colnam2);
%
tt = [t1 t2];
%
writetable(tt,xlsnam,'WriteVariableNames',false,'WriteMode','append');
%
return