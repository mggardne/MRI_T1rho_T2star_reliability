% 'AcquisitionDuration'
%
% Get MRI Series
%
sdir = dir('003');
sdirs = {sdir.name}';
sdirs = sdirs([sdir.isdir]');
sdirs = sdirs(startsWith(sdirs,'s'));
nsdir = size(sdirs,1);
%
%
%
info = dicominfo('003\s0501\i00016.dcm');
nams = fieldnames(info);
idt = find(contains(nams,'Time','IgnoreCase',true));
nams(idt)
idd = find(contains(nams,'Duration','IgnoreCase',true));
nams(idd)
ida = find(contains(nams,'Acq','IgnoreCase',true));
nams(ida)
%
info.InstanceCreationTime              % '094626.155'
info.StudyTime                         % '150751'
info.SeriesTime                        % '155015.89000'
info.AcquisitionTime                   % '155018.92'
info.AcquisitionDuration               % 111.9000
info.ContentTime                       % '155018.92'
info.RepetitionTime                    % 50
info.EchoTime                          % 0.4250
info.PerformedProcedureStepStartTime   % '150751'
info.PerformedProcedureStepEndTime     % '150751'
info.IssueTimeOfImagingServiceRequest  % '150751.428'

