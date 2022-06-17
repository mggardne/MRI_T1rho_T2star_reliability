xlsnam = 'MRIseries.xlsx';
%
dlst = ['MRIR 01-BC\Visit1\dicom_lst2.mat'
        'MRIR 01-BC\Visit2\dicom_lst2.mat'
        'MRIR 02-ES\Visit1\dicom_lst2.mat'
        'MRIR 02-ES\Visit2\dicom_lst2.mat'
        'MRIR 05-SC\Visit1\dicom_lst2.mat'
        'MRIR 05-SC\Visit2\dicom_lst2.mat'
        'MRIR 06-MM\Visit1\dicom_lst2.mat'
        'MRIR 06-MM\Visit2\dicom_lst2.mat'
        'MRIR 07-BP\Visit1\dicom_lst2.mat'
        'MRIR 07-BP\Visit2\dicom_lst2.mat'
        'MRIR 08-AS\Visit1\dicom_lst2.mat'
        'MRIR 08-AS\Visit2\dicom_lst2.mat'
        'MRIR 09-CR\Visit1\dicom_lst2.mat'
        'MRIR 09-CR\Visit2\dicom_lst2.mat'
        'MRIR 10-JF\Visit1\dicom_lst2.mat'
        'MRIR 10-JF\Visit2\dicom_lst2.mat'
        'MRIR 12-AO\Visit1\dicom_lst2.mat'
        'MRIR 12-AO\Visit2\dicom_lst2.mat'
        'MRIR 13-AK\Visit1\dicom_lst2.mat'
        'MRIR 13-AK\Visit2\dicom_lst2.mat'];
%
for ks = 1:20;
%
   load(dlst(ks,:),'etn','nimages','sn','stxt');
%
% Use Series Descriptions to Find T1rho Series
%
   idr = contains(stxt,'rho','IgnoreCase',true)&nimages>99;     % T1rho
   nr = sum(idr);
%
   stxtr = stxt(idr);        % T1rho series
%
% Use Series Descriptions to Find T2* Series
%
   ids = contains(stxt,'*')|contains(stxt,'T2','IgnoreCase',true);   % T2*
%
   ids = [ids; false];
   idds = diff(ids);
   idss = find(idds==1)+1;
   idse = find(idds==-1);
%
   ns = length(idss);        % Number of sets of T2* sequences
   nse = length(idse);
   if ns~=nse
     error([' *** ERROR in rd_dicom:  Number of starting and ending', ...
            ' indices for sets of T2* series are not equal!']);
   end
%
   ids = [idss idse];   % Index to start of sets (1st column) and end (2nd column)
%
   netn = diff(ids,1,2)+1;   % Number of echo times
%
clear idds idse idss nse;
%
% Check for Duplicate Echo Times
%
   netmin = min(netn);
%
   if netmin==max(netn)
     idd = 0;
   else
     idd = find(netn>netmin);          % Index to duplicates
   end
%
% Get Series Index and Echo Times as Matrices
%
   idse = zeros(netmin,ns);            % Series index matrix
%
   for k = 1:ns
%
      idsx = (ids(k,1):ids(k,2))';     % Index to series
      e = cell2mat(etn(idsx));         % Echo times for series
%
      if any(k==idd)                   % Correct duplicates
        [e,ide] = unique(e,'last');    % Get last duplicates
        if size(e,1)~=netmin
          error(' *** rd_dicom:  Number of echo times are not the same!');
        end
        idsx = idsx(ide);  % Correct series index
      end
%
      [~,idsrt] = sort(e);
      if all(idsrt~=(1:netmin)')
        error(' *** rd_dicom:  Echo times are not in order!');
      end
%
      if any(k==idd)                   % Warn about duplicates
        fprintf(1,'\n');
        warning(' *** WARNING:  Duplicate echo times!');
        fprintf(1,' Using series:\n');
        sstr = cellstr([int2str(sn(idsx)) repmat(' - ',netmin,1) ...
                        char(stxt{idsx})]);
        fprintf(1,'   Series %s\n',sstr{:})
        fprintf(1,'\n');
      end
%
      idse(:,k) = idsx;      % Series index into matrix
%
   end
   stxts = stxt(idse);       % T2* series
   stxts = stxts(1,:)';
%
   t = table(repmat(dlst(ks,:),4,1),stxtr,stxts);
   if ks==1
     writetable(t,xlsnam);
   else
     writetable(t,xlsnam,'WriteVariableNames',false, ...
                'WriteMode','append');
   end
%
end
%
return