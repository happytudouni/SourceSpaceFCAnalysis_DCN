function [] = BaselineDataProcessing_pps(participantnumber,age)
% The preprocessing program used for Xie et al.(2019),BMC Med; Xie, Toll, & Nelson (submitted to DCN);
% WX removed some personal information (e.g., directories in the personal and lab computers);
% WX also removed the scripts that were unrelated to preprocessing, e.g., loading files, creating event-list, etc.;
% because they were exclusively written for the Net Station data (in matlab format).
% Please contact WX if you are interested in those parts
% (wanze.xie@pku.edu.cn)
%-------------------------------------------------------------------------------------------------
plotarg    = 0;
replacearg = 0;
typeofdata = 0;
epochlength = 2; % define how long each epoch is (2 means 2s, could be changed to other numbers);

%% check if the final data already exists;
SegmentAverageFiles = '..\SegmentAverageFiles\';
if epochlength>1
Finaldataset = ['...' num2str(epochlength) 's.set'];
else
Finaldataset = ['...set'];
end
if exist([SegmentAverageFiles Finaldataset]) & replacearg == 0;
    disp([Finaldataset ' already exists! change replacearg to 1 if needs to replace it.']);
    return
end
%if ~exist('ALLEEG');eeglab;close;end
% eeglab; close; ft_defaults; run eeglab and ft_defaults only if needed (they haven't be run since the open of MATLAB);
global EEG;
addpath ..\programs\;
%% 1) load NetStation data


%% 2) Create Eventlist;


%% Remove the four channels with no data (only for HGSN for young children)
EEG = pop_select( EEG,'nochannel',{'E125' 'E126' 'E127' 'E128'});
EEG.urchanlocs(125:128) = [];

%% Filter
% EEGLAB FIR filter
EEG = pop_eegfiltnew(EEG,1,50);
% ERPLAB IIR filter
% EEG  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff', [1 50], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  8 );
if plotarg
pop_eegplot( EEG, 1, 1, 1);
end

%% Save the data as the EEGLAB format;
suffix = ['.' num2str(1) 'hzHighpass.set'];
tempname = strrep(filename,'.mat',suffix);
EEG  = pop_saveset( EEG, 'filename',tempname,'filepath', datapath); 


%% Find out extraordinarily bad channels using Faster;
list_properties = channel_properties(EEG, 1:EEG.nbchan, EEG.nbchan);
FASTbadIdx=min_z(list_properties); 
FASTbadChans=find(FASTbadIdx==1);

%% Segmentation;
binlistname = 'Binlist.txt';
EEG  = pop_binlister( EEG , 'BDF', ['...\' binlistname], 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG' ); 
EEG = eeg_checkset( EEG );
if plotarg
    pop_eegplot( EEG, 1, 1, 1);
end
segdura = 1000*segduration;
EEG = pop_epochbin( EEG , [0  segdura],  'none'); 
EEG = eeg_checkset( EEG );

%% Before running ICA, find out very bad channels with large amplitudes for more than 20% of the trials;
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

%% Run ICA
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'stop', 1E-7, 'interupt','off');

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


%% Channel interpolation in EEGLAB for epochs;
% Find artifaceted epochs by detecting outlier voltage after rejection of bad channels
EEG = pop_eegthresh(EEG,1, 1:numChans, -100, 100, EEG.xmin, EEG.xmax,0,0);
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
fieldtripfolder = '..\';
fieldtripname   = ['...mat'];

save([fieldtripfolder fieldtripname],'EEG_FT');

disp(['**************************** Done for ' num2str(participantnumber)  ' ****************************']);
