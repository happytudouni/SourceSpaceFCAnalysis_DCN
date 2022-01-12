function [] = BaselineDataProcessing(participantnumber,age)
% WX wrote this program in October, 2017; WX modified this on 6/12/2018 for KV's thesis project;
% WX modified this program in the first week of Nov, 2017 to make it also work for baseline data;
%-------------------------------------------------------------------------------------------------
plotarg    = 0;
replacearg = 0;
typeofdata = 0;
hz = 1;
epochlength = 1; % define how long each epoch is (1 means 1s, could be changed to other numbers);

%% check if the final data already exists;
SegmentAverageFiles = 'C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\data\LCN\BangladeshProject\SegmentAverageFiles\';
if epochlength>1
Finaldataset = ['Experiment 2 Subject ' num2str(participantnumber) ' Age ' num2str(age) ' 1Hz EEGLAB ArtifactsRejected_ICA_Edited_baseline' '_' num2str(epochlength) 's.set'];
else
Finaldataset = ['Experiment 2 Subject ' num2str(participantnumber) ' Age ' num2str(age) ' 1Hz EEGLAB ArtifactsRejected_ICA_Edited_baseline.set'];
end
if exist([SegmentAverageFiles Finaldataset]) & replacearg == 0;
    disp([Finaldataset ' already exists! change replacearg to 1 if needs to replace it.']);
    return
end
%if ~exist('ALLEEG');eeglab;close;end
% eeglab; close; ft_defaults; run eeglab and ft_defaults only if needed (they haven't be run since the open of MATLAB);
global EEG;
addpath C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\programs\;
%% 1) read the NSR filtered (0.3-30hz), segmented, baseline corrected Matlab .mat file
datapath = 'E:\data\BangladeshProject\NSRRawData\';
switch typeofdata
    case 'erp'
        typeofdata = 1;
        disp('This current process is on the ERP data!')
    case 'baseline'
        typeofdata = 0;
        disp('This current process is on the baseline data!')
end

if participantnumber<10
    participantnumberstring = ['00' num2str(participantnumber)];
elseif participantnumber<100
    participantnumberstring = ['0' num2str(participantnumber)];
else
    participantnumberstring = num2str(participantnumber);
end
if age == 24
    keystring = {participantnumberstring, '_2y.ns.mat'};
elseif age == 60
    keystring = {participantnumberstring, '_5y.ns.mat'};
else
    keystring = {participantnumberstring, 's.mat'};
end
filename = find_filename(keystring, datapath)
% Load in the raw matlab file;
fieldname = 'EEG_Segment1';
EEG = pop_importegimat([datapath filename],500.00, 0, fieldname); %This would remove the 129 channel that has 0s in its data;
% Re-calculate the 129 reference channel;
EEG.nbchan = 129;
%datach129 = zeros(1,size(EEG.data,2));
datach129 = mean(EEG.data,1);
EEG.data = [EEG.data;datach129];
EEG.epoch = 1;
EEG = readegilocs(EEG); % re-load the channel location for channel 129;
if plotarg
pop_eegplot( EEG, 1, 1, 1);
end


%% 2) Eventlist;
%The continuous data does not have an eventlist, but the epoch function in
%ERPLAB needs an eventlist to work, so create an eventlist for EEGLAB/ERPLAB;
[nch nsamples] = size(EEG.data);
nseconds = nsamples/500;
segduration = epochlength; %1 or 2s
nepochs  = nseconds/segduration;
EEG.epoch = 1;
item = 0;bepoch = 0;diff=0; dura=1000*segduration; enable=1;ecode = 1; label = 'TARG';
if segduration > 1
    erplistname = strrep(filename,'.mat',['.baseline.erplist_' num2str(epochlength) 's.txt']);
else
    erplistname = strrep(filename,'.mat','.baseline.erplist.txt');
end
if exist([datapath erplistname])
    delete([datapath erplistname])
end
item = 0;
for i = 1:nepochs;
    EEG.event(i).type     = 'target';
    EEG.event(i).latency  = 1+(i-1)*500*segduration; %latency in samples
    EEG.event(i).urevent  = 'target';
    EEG.event(i).duration = segduration;
    EEG.event(i).codes    = 1;
    onset = segduration*(i-1);
    item  = item +1;
    %string=sprintf('%s %u %s %u %u',EEG.event(i).type,EEG.event(i).latency,EEG.event(i).urevent,EEG.event(i).duration,EEG.event(i).codes);
    %disp(string);
    string=sprintf('%i %i %i %s %f %.2f %.1f 00000000 00000000 %i [    ]',item,bepoch,ecode,label,onset,diff,dura,enable);
    disp(string);
    dlmwrite([datapath erplistname],string,'-append','delimiter','');
end
    % Both EEGLAB and ERPLAB can do segmentation, but the artifacts detection functions are based on the epochs defined in ERLAB;
    % Also, for these functions to work, the filed of "EVENTLIST" needs to exist;


%% Remove the four channels and do average reference
EEG = pop_select( EEG,'nochannel',{'E125' 'E126' 'E127' 'E128'});
EEG.urchanlocs(125:128) = [];

%% Filter
EEG  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff', [1 50], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  8 );
%EEG1  = pop_saveset( EEG1, 'filename',tempname,'filepath', datapath); %baseline data have not been filetered yet;
if plotarg
pop_eegplot( EEG, 1, 1, 1);
end

%% Save the data as the EEGLAB format;
suffix = ['.' num2str(1) 'hzHighpass.set'];
tempname = strrep(filename,'.mat',suffix);
EEG  = pop_saveset( EEG, 'filename',tempname,'filepath', datapath); %baseline data have not been filetered yet;

%% Re-load the data;
EEG = pop_loadset('filename',tempname,'filepath',datapath);
chanlocs = EEG.chanlocs;
nbchan = EEG.nbchan;
if nbchan == 125;
    nbchan = nbchan -1;
    chanlocs(end) = [];
end
% import eventlist;
if segduration>1
    keystring = {num2str(participantnumber),['baseline.erplist_' num2str(epochlength) 's.txt']};
else
    keystring = {num2str(participantnumber),'baseline.erplist.txt'};
end
elistname = find_filename(keystring, datapath);
EEG = pop_importeegeventlist( EEG, [datapath elistname], 'ReplaceEventList', 'on' );
EEG = eeg_checkset( EEG);

%% Find out extraordinarily bad channels using Faster;
list_properties = channel_properties(EEG, 1:EEG.nbchan, EEG.nbchan);
FASTbadIdx=min_z(list_properties); % WX changed the threshold from 3 std to 4 std in FASTER;
FASTbadChans=find(FASTbadIdx==1);

%% Import Binlist and Do Segmentation;
binlistname = 'BangladeshBaselineBinList.txt';
EEG  = pop_binlister( EEG , 'BDF', ['C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\programs\BangladeshProject\' binlistname], 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG' ); % GUI: 24-Oct-2017 15:26:01
EEG = eeg_checkset( EEG );
if plotarg
pop_eegplot( EEG, 1, 1, 1);
end
segdura = 1000*segduration;
EEG = pop_epochbin( EEG , [0  segdura],  'none'); % GUI: 24-Oct-2017 15:27:00
EEG = eeg_checkset( EEG );

%% Find out additionally bad channels with large amplitudes for more than 20% of the trials;
vol_thrs = [-500 500]; % [lower upper] threshold limit(s) in mV.
emg_thrs = [-100 30]; % [lower upper] threshold limit(s) in dB.
emg_freqs_limit = [20 50]; % [lower upper] frequency limit(s) in Hz.

numChans =EEG.nbchan; % Find the number of channels
numEpochs =EEG.trials; % Find the number of epochs
chanCounter = 1;
badChans = [];
badEpochOutput = [];
for ch=1:numChans
    % Find artifaceted epochs by detecting outlier voltage
    EEG = pop_eegthresh(EEG,1, ch, vol_thrs(1), vol_thrs(2), EEG.xmin, EEG.xmax, 0, 0);
    EEG = eeg_checkset( EEG );

    % 1         : data type (1: electrode, 0: component)
    % 0         : Display with previously marked rejections? (0: no, 1: yes)
    % 0         : Reject marked trials? (0: no, 1:yes)

    % Find artifaceted epochs by using thresholding of frequencies in
    % the data. This method mainly rejects muscle movement (EMG) artifacts
    EEG = pop_rejspec( EEG, 1,'elecrange',ch ,'method','fft','threshold', emg_thrs, 'freqlimits', emg_freqs_limit, 'eegplotplotallrej', 0, 'eegplotreject', 0);

    % method                : method to compute spectrum (fft)
    % threshold             : [lower upper] threshold limit(s) in dB.
    % freqlimits            : [lower upper] frequency limit(s) in Hz.
    % eegplotplotallrej     : 0 = Do not superpose rejection marks on previous marks stored in the dataset.
    % eegplotreject         : 0 = Do not reject marked trials (but store the  marks.

    % Find number of artifacted epochs
    EEG = eeg_checkset( EEG );
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    artifacted_epochs=EEG.reject.rejglobal;

    % Find bad channel / channel with more than 33.3% artifacted epochs
    if sum(artifacted_epochs) > (numEpochs*0.333)
        badChans(chanCounter) = ch;
        badEpochOutput(chanCounter) = sum(artifacted_epochs);
        chanCounter=chanCounter+1;
    end
end

%% Save the bad channel information;
badChans = union(FASTbadChans,badChans);
badChans(badChans == 125) = []; %the Cz (channel 125) will be removed anyway later;
if segduration>1
    subject_name=['Participantnumber ' num2str(participantnumber) ' Age ' num2str(age) '_Bad_Channels_' num2str(epochlength) 's.mat'];
else
    subject_name=['Participantnumber ' num2str(participantnumber) ' Age ' num2str(age) '_Bad_Channels.mat'];
end
save([datapath subject_name],'badChans');
if length(badChans) > 18;
    disp(['Participant ' num2str(participantnumber) ' has more than 18 bad channels, so should not be used.']);
    return
end
%% Remove the bad channels and Cz;
chan2Bremoved = union(badChans,125);
EEG = pop_select( EEG,'nochannel',chan2Bremoved );

%% Reject very bad epochs;
% Find the artifacted epochs across all channels and reject them before doing ICA.
% Find the number of channels again since some channels have been removed
numChans =EEG.nbchan;

% Find artifaceted epochs by detecting outlier voltage after rejection of bad channels
EEG = pop_eegthresh(EEG,1, 1:numChans, vol_thrs(1), vol_thrs(2), EEG.xmin, EEG.xmax,0,0);
EEG = eeg_checkset( EEG );

% Find artifaceted epochs by using thresholding of frequencies in the data.
EEG = pop_rejspec( EEG, 1,'elecrange', 1:numChans , 'method', 'fft', 'threshold', emg_thrs ,'freqlimits', emg_freqs_limit, 'eegplotplotallrej', 0, 'eegplotreject', 0);
%[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Find the number of artifacted epochs and reject them
EEG = eeg_checkset( EEG );
EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
% if there are more than 5 bad channels then reject the epoch; otherwise, interpolate the bad channels;
badchannels_firstround = EEG.reject.rejglobalE; % save this information for artifacts detection and channel interpolation down below;
badchs4epochs= sum(EEG.reject.rejglobalE,1);
reject_artifacted_epochs=find(badchs4epochs>5);
badchannels_firstround(:,reject_artifacted_epochs) = [];
badchmatrix = EEG.reject.rejglobalE;
% interpolate for channels; Loop through each epoch, select it, run interp, save data
numEpochs = EEG.trials;
tmpData = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
for e = 1:numEpochs
    % Initialize variables EEGe and EEGe_interp;
    if ismember(e,reject_artifacted_epochs) % this epoch has more than 5 bad channels;
    else
        EEGe = [];
        EEGe_interp = [];
        badChanNum = [];

        %% select only this epoch (e)
        EEGe = pop_selectevent( EEG, 'epoch', e, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        badChanNum = find(badchmatrix(:,e)); % find which channels are bad for this epoch
        EEGe_interp = eeg_interp(EEGe, badChanNum); % interpolate the bad chans for this epoch
        tmpData(:,:,e) = EEGe_interp.data; % store interpolated data into matrix
        EEG.reject.rejglobal(e)=0;
        EEG.reject.rejglobalE(badChanNum,e) = 0;
    end
end
EEG.data = tmpData;
% reject the bad epochs;
EEG = pop_rejepoch( EEG, reject_artifacted_epochs ,0);
% save rejected epochs;
if segduration>1
    rejectedepochs = ['Participantnumber ' num2str(participantnumber) ' Age ' num2str(age) '_RejectedEpochs_' num2str(epochlength) 's'];
else
    rejectedepochs = ['Participantnumber ' num2str(participantnumber) ' Age ' num2str(age) '_RejectedEpochs'];
end
save([datapath rejectedepochs],'reject_artifacted_epochs');
% Give a name to the copied dataset and save after cleaning
Copied_Data=['Participantnumber ' num2str(participantnumber) ' Age ' num2str(age) '_ready_for_ICA_Run'];
%EEG = pop_saveset( EEG, 'filename',Copied_Data, 'filepath', datapath);

%% Run ICA
if segduration>1
    Copied_ICA_Data = ['Participantnumber ' num2str(participantnumber) ' Age ' num2str(age) '_ICA_Run_' num2str(epochlength) 's.set'];
else
    Copied_ICA_Data = ['Participantnumber ' num2str(participantnumber) ' Age ' num2str(age) '_ICA_Run.set'];
end
if exist([datapath Copied_ICA_Data]) & replacearg == 0;
%    EEG = pop_loadset('filename',Copied_ICA_Data,'filepath',datapath);
else
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'stop', 1E-7, 'interupt','off');
    % Save the ICA copied data and ICA weights;
    EEG = pop_saveset( EEG, 'filename',Copied_ICA_Data,'filepath', datapath);
end

%% find and save ICA weights;
icainfo.ICA_WINV=EEG.icawinv;
icainfo.ICA_SPHERE=EEG.icasphere;
icainfo.ICA_WEIGHTS=EEG.icaweights;
icainfo.ICA_CHANSIND=EEG.icachansind;
if segduration>1
    icainfoname = ['Participantnumber ' num2str(participantnumber) ' Age ' num2str(age) '_ICA_Weights_' num2str(epochlength) 's.mat'];
else
    icainfoname = ['Participantnumber ' num2str(participantnumber) ' Age ' num2str(age) '_ICA_Weights.mat'];
end
save([datapath icainfoname],'icainfo');

%% calculate the two EOG channels;
% replace the bad channels in case some of them are EOG channels;
% however, these channels are not replaced yet for the curreng EEG;
nbchan = 124; %This is the channel number with all electrodes;
tempdata = zeros(nbchan,EEG.pnts,EEG.trials,'single');
goodChans = ~ismember([1:nbchan],badChans);
tempdata(goodChans,:,:) = EEG.data;
EEG1 = EEG;
EEG1.data = tempdata;
EEG1.chanlocs = chanlocs;
EEG1.nbchan = nbchan;
EEG1 = eeg_interp(EEG1,badChans);
VEOG = mean(EEG1.data([14 15 21],:,:),1);
HEOG = mean(EEG1.data([25 32 26],:,:),1) - mean(EEG1.data([1 2 8],:,:),1);

EEG.data(end+1,:,:) = VEOG;
EEG.data(end+1,:,:) = HEOG;
EEG.nbchan = size(EEG.data,1);
EEG.chanlocs(end+1).labels = 'VEOG';
EEG.chanlocs(end+1).labels = 'HEOG';
EEG = eeg_checkset( EEG );

%% Define the badICs with my grogram calling SASICA and Adjust;
cfg = [];cfg.plotarg = 0;
[badICs EEG] = DefineArtificialICAs_AdjustAndSASICA(cfg,EEG); % if include EEG as an output structure then it will be the one after bad ICs removed;
% Remove the two EOG channels;
EEG = pop_select( EEG,'nochannel',{'VEOG' 'HEOG'});

%% Reject the bad components
%EEG.reject.gcompreject(badICs) = 1;
%EEG = eeg_checkset( EEG );
%EEG = pop_subcomp( EEG, badICs, 0);

%% Channel interpolation in EEGLAB for epochs;
% Find artifaceted epochs by detecting outlier voltage after rejection of bad channels
if segduration>1
    if age < 60
        EEG = pop_eegthresh(EEG,1, 1:numChans, -150, 150, EEG.xmin, EEG.xmax,0,0);
    else
        EEG = pop_eegthresh(EEG,1, 1:numChans, -150, 150, EEG.xmin, EEG.xmax,0,0); % later, may think about a way to make the threshold more rigorous;
    end
else
    EEG = pop_eegthresh(EEG,1, 1:numChans, -100, 100, EEG.xmin, EEG.xmax,0,0);
end
EEG = eeg_checkset( EEG );

% Find artifaceted epochs by using thresholding of frequencies in the data.
EEG = pop_rejspec( EEG, 1,'elecrange', 1:numChans , 'method', 'fft', 'threshold', emg_thrs ,'freqlimits', emg_freqs_limit, 'eegplotplotallrej', 0, 'eegplotreject', 0);

%% Find the number of artifacted epochs and reject them
EEG = eeg_checkset( EEG );
EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
% if there are more than 5 bad channels then reject the epoch; otherwise, interpolate the bad channels;
% include the first round bad channels;
badchannels_secondround = EEG.reject.rejglobalE;
banchannels_1st2ndrounds = badchannels_firstround + badchannels_secondround;
banchannels_1st2ndrounds(banchannels_1st2ndrounds == 2) = 1;
badchs4epochs= sum(banchannels_1st2ndrounds,1);
%badchs4epochs= sum(EEG.reject.rejglobalE,1);
maxchan2breplaced = 18 - length(badChans);
reject_artifacted_epochs=find(badchs4epochs>maxchan2breplaced);
badchmatrix = EEG.reject.rejglobalE;
% interpolate for channels; Loop through each epoch, select it, run interp, save data
numEpochs = EEG.trials;
tmpData = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
for e = 1:numEpochs
    % Initialize variables EEGe and EEGe_interp;
    if ismember(e,reject_artifacted_epochs) % this epoch has more than 5 bad channels;
    else
        EEGe = [];
        EEGe_interp = [];
        badChanNum = [];

        %% select only this epoch (e)
        EEGe = pop_selectevent( EEG, 'epoch', e, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        badChanNum = find(badchmatrix(:,e)); % find which channels are bad for this epoch
        EEGe_interp = eeg_interp(EEGe, badChanNum); % interpolate the bad chans for this epoch
        tmpData(:,:,e) = EEGe_interp.data; % store interpolated data into matrix
        EEG.reject.rejglobal(e)=0;
        EEG.reject.rejglobalE(badChanNum,e) = 0;
    end
end
EEG.data = tmpData;
% reject the bad epochs;
EEG = pop_rejepoch( EEG, reject_artifacted_epochs ,0);

%% Interpolate the overall bad channels that were removed from the beginning;
tempdata = zeros(nbchan,EEG.pnts,EEG.trials,'single');
goodChans = ~ismember([1:nbchan],badChans);
tempdata(goodChans,:,:) = EEG.data;
EEG.data = tempdata;
EEG.chanlocs = chanlocs;
EEG.nbchan = nbchan;
EEG = eeg_interp(EEG,badChans);

%% Re-referencing to average;
EEG = pop_reref(EEG, []);

%% Save the data in EEGLAB and Fieldtrip format;
EEG  = pop_saveset(EEG, 'filename',Finaldataset,'filepath', SegmentAverageFiles);

EEG_FT = eeglab2fieldtrip(EEG,'preprocessing','none');
fieldtripfolder = 'C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\data\LCN\BangladeshProject\Fieldtrip\rawdata\';
if segduration > 1
    fieldtripname = ['Experiment 2 Subject ' num2str(participantnumber) ' Age ' num2str(age) ' Fieldtrip 1Hz Highpass_New_' num2str(epochlength) 's.mat'];
else
    fieldtripname = ['Experiment 2 Subject ' num2str(participantnumber) ' Age ' num2str(age) ' Fieldtrip 1Hz Highpass_New.mat'];
end
save([fieldtripfolder fieldtripname],'EEG_FT');

disp(['**************************** Done for ' num2str(participantnumber)  ' ****************************']);
