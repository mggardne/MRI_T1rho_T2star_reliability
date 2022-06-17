%#######################################################################
%
%                     * DICOM HeaDeR Fix Program *
%
%          M-File which reads the original T2* DICOM images from a
%     series and the resized T2* DICOM images.  The range of image data
%     is compared and output to the screen.  The scaling factors from
%     the original series are written to the resized series.
%
%     NOTES:  1.  For use with subject 08-AS on Visit 2 for the T2* for
%             the unloaded left leg.
%
%     04-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Directories with T2* DICOM Files
%
ddirs = [ 's0401';      % Original series
          's3801' ];    % Resized series
%
% DICOM File Names
%
fnams = [ 'i00001.dcm'
          'i00002.dcm'
          'i00003.dcm'
          'i00004.dcm'
          'i00005.dcm'
          'i00006.dcm'
          'i00007.dcm'
          'i00008.dcm'
          'i00009.dcm'
          'i00010.dcm'
          'i00011.dcm'
          'i00012.dcm'
          'i00013.dcm'
          'i00014.dcm'
          'i00015.dcm'
          'i00016.dcm'
          'i00017.dcm'
          'i00018.dcm'
          'i00019.dcm'
          'i00020.dcm'
          'i00021.dcm'
          'i00022.dcm'
          'i00023.dcm'
          'i00024.dcm'
          'i00025.dcm'
          'i00026.dcm'
          'i00027.dcm'
          'i00028.dcm'
          'i00029.dcm'
          'i00030.dcm'
          'i00031.dcm'
          'i00032.dcm'
          'i00033.dcm'
          'i00034.dcm'
          'i00035.dcm'
          'i00036.dcm'
          'i00037.dcm'
          'i00038.dcm'
          'i00039.dcm'
          'i00040.dcm'
          'i00041.dcm'
          'i00042.dcm'
          'i00043.dcm'
          'i00044.dcm'
          'i00045.dcm'
          'i00046.dcm'
          'i00047.dcm'
          'i00048.dcm' ];
%
% Loop through the Files
%
frmt = [' %2i   %5.2f      %1i         %1i      %4i        %1i ', ...
        '       %4i       %3i\n'];
%
hdr = {'File'; 'Slope'; 'Intercept'; 'OrigMin'; 'OrigMax'; ...
       'ResizeMin'; 'ResizeMax'; 'DiffMax'};
fprintf(1,'\n%s',hdr{1});
for k = 2:7
   fprintf(1,'  %s',hdr{k});
end
fprintf(1,'  %s\n',hdr{8});
%
for k = 1:48
   info = dicominfo(fullfile(ddirs(1,:),fnams(k,:)));
   img = dicomread(info);
   mn1 = double(min(img(:)));
   mx1 = double(max(img(:)));
   sl = double(info.RescaleSlope);
   b = double(info.RescaleIntercept);
%
   info = dicominfo(fullfile(ddirs(2,:),fnams(k,:)));
   img = dicomread(info);
   mn2 = double(info.SmallestImagePixelValue);
   mx2 = double(info.LargestImagePixelValue);
%
   fprintf(1,frmt,[k sl b mn1 mx1 mn2 mx2 mx1-mx2]);
%
   info.RescaleSlope = single(sl);
   info.RescaleIntercept = single(b);
   info.RescaleType = 'normalized';
%
   dicomwrite(img,fullfile(ddirs(2,:),fnams(k,:)),info, ...
              'CreateMode','copy');
%
end