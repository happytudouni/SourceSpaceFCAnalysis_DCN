%% 0. Preamble.

% This script describes the preprocessing procedures to reduce artifacts and prepare the EEG for
% connectivity analyses. This script is annotated with comments throughout to aid in clarity.

% Conventions:
% functions prefixed by RTT_ are custom functions
% variables prefixed by thresh are values that are manually set specifically for these analyses
% fields in the EEG structure beginning with a capital letter are custom fields (ex. EEG.ChansBad is a custom field and
% EEG.data is a standard field)

clearvars
close all
clc

%% 1. Convert acquisitions to .set format

% This stage simply converts the raw acquisition data into the .set format used by EEGLAB. An empty EEGLAB basic EEG
% structure is loaded, populated with the corresponding data for each field, and then saved as a .set file to the
% Processing Stage 0 folder (PROC0).

% folder containing acquisition data
folderAcquisition = 'X:\COLLAB\BABY\EEG\RAW';
% folder to output preprocessed data
folderPreProc0 = 'R:\COLLAB\BABY\EEG\PREPROC\PREPROC0';
% RTT_DirList returns a table of the directory's file listing
D = RTT_DirList(folderAcquisition);
% restrict table to the acquisition files (.mat)
D = D(endsWith(D.PathAndName,'.mat'),:);
% generate uniform naming convention
D.Subject = cellfun(@char,regexp(D.Name,'Subject (\d+)','tokens','once'),'UniformOutput',false);
D.Subject = strcat('00',D.Subject);
D.Subject = cellfun(@(x) x(end-2:end),D.Subject,'UniformOutput',false);
D.Age = cellfun(@char,regexp(D.Name,'Age (\d+)','tokens','once'),'UniformOutput',false);
D.NameNew = strcat('A',D.Age,'_S',D.Subject,'.set');
% generate full path to the .set file in the Processing Stage 0 folder (PROC0)
D.PathAndNameNew = fullfile(folderPreProc0,D.NameNew);
for iFile=1:height(D)
   % check to see if this file has already been processed to prevent overwriting when processing newly added data
   if ~exist(D.PathAndNameNew{iFile},'file')
      % load eeglab empty shell
      load('R:\COLLAB\BABY\EEG\ADMIN\masterEmptyGlobalsEEGLAB.mat')
      % load fieldtrip acquisition data
      load(D.PathAndName{iFile});
      % assign fieldtrip acquisition data into EEG structure
      data = horzcat(EEG_FT.trial{:});
      times = horzcat(EEG_FT.time{:});
      srate = EEG_FT.fsample;
      chanlocs = EEG_FT.chanlocs;
      EEG.subject = extractBefore(D.NameNew{iFile},'.set');
      EEG.data = data;
      EEG.times = times;
      EEG.srate = srate;
      EEG.chanlocs = chanlocs;
      [EEG,changes] = eeg_checkset(EEG);
      % save .set file
      EEG = pop_saveset(EEG,'filename',D.PathAndNameNew{iFile},'savemode','onefile','version','7.3');
   end
end
disp('section 1 complete')

%% 2. Assign uniform channel structure compatible with Brainstorm 

clearvars
clc

% standard channel organization structure derived from Brainstorm fitted cap templates
load('R:\COLLAB\BABY\EEG\ADMIN\masterChaninfo.mat')
load('R:\COLLAB\BABY\EEG\ADMIN\masterChanlocsTable.mat')

% declare the type of cap used in these acquisitions
threshCap = 'EGI128'; % 128 channel EGI saline geodesic cap (124 channels used)

folderPreProc0 = 'R:\COLLAB\BABY\EEG\PREPROC\PREPROC0';
D = RTT_DirList(folderPreProc0);
D.PathAndNameNew = strrep(D.PathAndName,'PROC0','PROC1');

for iFile=1:height(D)
   % check to see if this file has already been processed to prevent overwriting when processing newly added data
   if ~exist(D.PathAndNameNew{iFile},'file')
      EEG = pop_loadset(D.PathAndName{iFile});
      chanlocs = struct2table(EEG.chanlocs);
      chanlocs.labels = upper(chanlocs.labels);
      % generate uniform channel labels
      chanNumChar = regexp(chanlocs.labels,'\d+','match','once');
      chanNumChar = strcat('000',chanNumChar);
      chanNumChar = cellfun(@(x) x(end-2:end),chanNumChar,'UniformOutput',false);
      chanlocs.labels = strcat('E',chanNumChar);
      % obtain uniform channel structure from joining to master table
      chanlocs.IndexOriginal = (1:height(chanlocs))';
      chanlocs = chanlocs(:,{'IndexOriginal','labels'});
      chanlocs = outerjoin(chanlocs,masterChanlocsTable(ismember(masterChanlocsTable.Cap,threshCap),:),'MergeKeys',true,'Type','left');
      % reorder data in accordance with uniform channel structure
      chanlocs = chanlocs(~isnan(chanlocs.IndexChannel),:);
      chanlocs = sortrows(chanlocs,'IndexChannel');
      chanlocs.IndexChannel = (1:height(chanlocs))';
      EEG.data = EEG.data(chanlocs.IndexOriginal,:);
      % update chanlocs in the EEG structure
      chanlocs = movevars(chanlocs,{'IndexOriginal','IndexCap','Cap'},'After','ColorB');
      chanlocs = movevars(chanlocs,'labels','After','IndexChannel');
      EEG.chanlocs = table2struct(chanlocs);
      EEG.chaninfo = masterChaninfo;
      EEG.nbchan = height(chanlocs);
      % save .set file
      EEG = pop_saveset(EEG,'filename',D.PathAndNameNew{iFile},'savemode','onefile','version','7.3');
   end
end
disp('section 2 complete')

%% 3. Calculate measures for quality inspection by a human expert

clearvars
clc

folderPreProc1 = 'R:\COLLAB\BABY\EEG\PREPROC\PREPROC1';
folderPreProc2 = strrep(folderPreProc1,'PROC1','PROC2');
D = RTT_DirList(folderPreProc1);
% .set file with notch filtering complete
D.PathAndNameNew = strrep(D.PathAndName,'\PROC1\','\PROC2\');
% inspection data for human review
D.PathAndNameNewInspect = strrep(D.PathAndNameNew,'.set','_INSPECT.mat');
% placeholder to prevent overwriting when running multiple instances of this script
D.PathAndNameNewPlaceholder = strrep(D.PathAndNameNewInspect,'_INSPECT.mat','_INPROGRESS.mat');
D.Subject = extractBefore(D.name,'.set');
threshDBnumIQRs = 20; % decibel threshhold to define a channel as clearly bad so it doesn't skew the scale during manual inspection

for iFile=1:height(D)
   % check to see if the files for this subject already exist
   if isempty(dir([folderPreProc2,'\',D.Subject{iFile},'*']))
      % save placeholder to prevent overwriting
      save(D.PathAndNameNewPlaceholder{iFile},'iFile')
      EEG = pop_loadset(D.PathAndName{iFile});
      % verify data is in double precision
      EEG.data = double(EEG.data);
      % calculate raw spectra
      [spectraRaw,freqs] = spectopo(EEG.data, EEG.pnts, EEG.srate,'winsize',EEG.srate,'overlap',EEG.srate/2,'percent',100,'plot','off');
      % define notch table
      notchTable = array2table((60:60:EEG.srate/2)','VariableNames',{'Peak'});
      notchTable.LoStop = notchTable.Peak - 4;
      notchTable.HiStop = notchTable.Peak + 4;
      % filter line noise harmonics
      for j=1:height(notchTable)
         EEG = pop_eegfiltnew(EEG,notchTable.LoStop(j),notchTable.HiStop(j),[],1);
      end
      % calculate notch spectra
      [spectraNotch,~] = spectopo(EEG.data, EEG.pnts, EEG.srate,'winsize',EEG.srate,'overlap',EEG.srate/2,'percent',100,'plot','off');
      EEGnotch = EEG;
      % calculate inspection data
      % note eeglab filtering functions take -6 dB (mid-transition band) as the input parameter
      % the desired passband cutoff is denoted to the right of the function
      % low pass
      EEG = pop_eegfiltnew(EEG,0,40); % cutoff 45 Hz
      % high pass
      EEG = pop_eegfiltnew(EEG,2,0); % cutoff 1 Hz
      % calculate spectra in physiological spectrum
      [spectraPhysio,~] = spectopo(EEG.data, EEG.pnts, EEG.srate,'winsize',EEG.srate,'overlap',EEG.srate/2,'percent',100,'plot','off');
      % determine initial bad chans by extreme spectra and record in chanlocs table
      freqIndex = freqs == 10;
      chansBad = find(abs((spectraPhysio(:,freqIndex) - median(spectraPhysio(:,freqIndex))) ./ iqr(spectraPhysio(:,freqIndex))) > threshDBnumIQRs);
      chanlocsTable = struct2table(EEG.chanlocs);
      for r=1:max(chanlocsTable.IndexRegion)
         chanlocsTable.IsSeed(find(chanlocsTable.IndexRegion == r,1,'first')) = true;
      end
      chanlocsTable.IsBad(chansBad) = true;
      chanlocsTable.ReasonBad = repmat({''},[height(chanlocsTable),1]);
      chanlocsTable.ReasonBad(chansBad) = {'spectra'};
      % resample temporary inspection time series to reduce storage requirement
      EEG = pop_resample(EEG,50);
      % data for inspection figure
      t = (1:EEG.pnts) ./ EEG.srate;
      data = single(EEG.data);
      % save inspection data
      save(D.PathAndNameNewInspect{iFile},'spectra*','freqs','*Table','t','data')
      % assign custom fields to EEGnotch .set
      EEG = EEGnotch;
      EEG.Thresh.threshDBnumIQRs = threshDBnumIQRs;
      EEG.NotchTable = notchTable;
      EEG.Spectra.Raw.Freqs = freqs;
      EEG.Spectra.Raw.Spectra = spectraRaw;
      EEG.Spectra.Notch.Freqs = freqs;
      EEG.Spectra.Notch.Spectra = spectraNotch;
      EEG.Spectra.Physio.Freqs = freqs;
      EEG.Spectra.Physio.Spectra = spectraPhysio;
      EEG.ChansBad = chansBad;
      EEG.chanlocs = table2struct(chanlocsTable);
      % save EEG with notch filtered data
      EEG = pop_saveset(EEG,'filename',D.PathAndNameNew{iFile},'savemode','onefile','version','7.3');
      clearvars -except D folder* thresh* iFile
      % delete placeholder
      delete(D.PathAndNameNewPlaceholder{iFile})
   end
end

disp('section 3 complete')

%% 4. Manually inspect the data to identify bad channels

% This stage presents a figure for inspection. The figure displays the time series at top in butterfly format by head
% region. The middle portion presents the raw, notched, and physiological power spectra. The bottom panels display the
% physiological power spectra by head region. The inspection cycle is as follows:
% 1. Display inspection figure. Figure title (middle, right) is red, indicating the script is ready for human input.
% 2. Left click on any spectra to identify as bad channel. If a mistake is made, right click anywhere to restart this
% iteration.
% 3. Left click on the time series when done. Figure title will turn yellow while calculations update. Wait until it
% returns to red to proceed.
% 4. Repeat steps 2-3 as necessary. When fully complete left click again on time series. Figure title will turn to 
% magenta. Left click to complete this subject, the figure title turns black, and the figure and inspection files are
% saved. Right click if there are issues needing to be addressed and the filename will be changed to !REINSPECT.
% 5. Cycle this section until all subjects processed. Manually closing the figure will cause an error and halt the
% script.

clearvars
clc

folderPreProc2 = 'R:\COLLAB\BABY\EEG\PREPROC\PREPROC2';
D = RTT_DirList(folderPreProc2);
D = D(endsWith(D.name,'_INSPECT.mat'),:);
% completed inspection from human review
D.PathAndNameNewInspected = strrep(D.PathAndName,'_INSPECT.mat','_INSPECTED.mat');
D.PathAndNameNewFigure = strrep(D.PathAndNameNewInspected,'.mat','.jpg');

for iFile=1:height(D)
   if ~exist(D.PathAndNameNewInspected{iFile},'file')
      load(D.PathAndName{iFile})

      nbchan = height(chanlocsTable);
      numRegions = max(chanlocsTable.IndexRegion);
      chanColors = [chanlocsTable.ColorR,chanlocsTable.ColorG,chanlocsTable.ColorB] ./ 255;
      chanColorsRegions = zeros(numRegions,3);
      chansBad = find(chanlocsTable.IsBad);

      data(chansBad,:) = NaN;
      spectraRaw(chansBad,:) = NaN;
      spectraNotch(chansBad,:) = NaN;
      spectraPhysio(chansBad,:) = NaN;

      scale = 10 * iqr(data(:));
      dataRB = data;
      for r=1:numRegions
         indices = chanlocsTable.IndexRegion == r;
         if any(indices)
            chanColorsRegions(r,:) = chanColors(indices & chanlocsTable.IsSeed,:);
         end
         dataRB(indices,:) = data(indices,:) + ((numRegions - r) * 2 * scale);
      end
      % fixed variables

      threshIndicesPhysio = freqs <= 40;

      threshMonitorSelected = 1;
      threshTimeSecondsMax = 0;
      threshFigRes = 72; % dpi
      threshTmax = 0;
      threshFigTitle = strrep(extractBefore(D.name{iFile},'_INSPECT.mat'),'_','\_');

      % bands
      threshBandDelta = [0,3.5];
      threshBandTheta = [3.5,7.5];
      threshBandAlpha = [7.5,12.5];
      threshBandBeta1 = [12.5,22.5];
      threshBandBeta2 = [22.5,30.5];
      threshBandGamma = [30.5,45];

      % inspection plot function
      recompute = true;
      while recompute
         recompute = false;
         % Initialize Figure
         if exist('fig','var') && isvalid(fig)
            clf(fig)
         else
            set(0,'Units','normalized')
            monitorPositions = get(0,'MonitorPositions');
            fig = figure;
            fig.Tag = 'fig';
            fig.MenuBar = 'none';
            fig.ToolBar = 'none';
            fig.Name = 'Proc 2 Bad Chans (Spectra) Visual Inspection of TS and Spectra';
            fig.Color = 'w';
            fig.Resize = 'off';
            fig.Units = 'normalized';
            if threshMonitorSelected == 0
               fig.Visible = 'off';
               fig.OuterPosition = monitorPositions(1,:);
            else
               fig.Visible = 'on';
               fig.OuterPosition = monitorPositions(threshMonitorSelected,:);
            end
         end

         % Text Output
         axTX = axes;
         axTX.Tag = 'axTX';
         axTX.Toolbar.Visible = 'off';
         axTX.Units = 'normalized';
         axTX.Box = 'off';
         axTX.Color = 'none';
         axTX.Position = [0,0,1,1];
         axTX.XLim = [0,1];
         axTX.YLim = [0,1];
         axTX.XAxis.Visible = 'off';
         axTX.YAxis.Visible = 'off';

         txSPVP = text;
         txSPVP.Parent = axTX;
         txSPVP.String = threshFigTitle;
         txSPVP.FontSize = 18;
         txSPVP.FontWeight = 'bold';
         txSPVP.HorizontalAlignment = 'right';
         txSPVP.VerticalAlignment = 'top';
         txSPVP.Position = [1,.6,0];
         txSPVP.Color = 'k';

         txBC = text;
         txBC.Parent = axTX;
         if isempty(chansBad)
            txBC.String = sprintf('Bad Chans: %u (%g%%)\n',length(chansBad),round(100 * length(chansBad) / nbchan,1));
         else
            txBC.String = sprintf('Bad Chans: %u (%g%%)\n%s',length(chansBad),round(100 * length(chansBad) / nbchan,1),strtrim(regexprep(mat2str(chansBad),'\D','\n')));
         end
         txBC.FontSize = 8;
         txBC.FontWeight = 'bold';
         txBC.HorizontalAlignment = 'right';
         txBC.VerticalAlignment = 'top';
         txBC.Position = [1,.5725,0];
         txBC.Color = 'r';

         % Regional Butterfly
         axRB = axes;
         axRB.Tag = 'axRB';
         axRB.Toolbar.Visible = 'off';
         axRB.Units = 'normalized';
         axRB.Box = 'off';
         axRB.Color = 'none';
         x1RB = .01;
         x2RB = .995;
         y1RB = .6;
         y2RB = 1;
         wRB = x2RB - x1RB;
         hRB = y2RB - y1RB;
         axRB.InnerPosition = [x1RB, y1RB, wRB, hRB];
         axRB.Box = 'off';
         for r=1:numRegions
            indices = find(chanlocsTable.IndexRegion == r);
            dataRBregion = dataRB(indices,:);
            [~,j] = sort(iqr(dataRBregion,2),'descend');
            dataRBregion = dataRBregion(j,:);
            hold on
            for j=1:size(dataRBregion,1)
               plot(axRB,t,dataRBregion(j,:))
            end   
         end
         axRB.XTick = [];
         axRB.YTick = [];
         if threshTmax == 0
            threshTmax = t(end);
         end
         axRB.XLim = [0,threshTmax];
         axRB.YLim = [-scale,2 * (numRegions - 1) * scale + scale];
         axRB.XAxis.Visible = 'off';
         axRB.YAxis.Visible = 'off';

         axRR = axes;
         axRR.Tag = 'axRR';
         axRR.Toolbar.Visible = 'off';
         axRR.Units = 'normalized';
         x1RR = 0;
         x2RR = x1RB;
         y1RR = y1RB;
         y2RR = y2RB;
         wRR = x2RR - x1RR;
         hRR = y2RR - y1RR;
         axRR.InnerPosition = [x1RR, y1RR, wRR, hRR];
         axRR.Box = 'off';
         axRR.XTick = [];
         axRR.YTick = [];
         axRR.XLim = [0,1];
         axRR.YLim = [0,1];
         for r=1:numRegions
            index = find(chanlocsTable.IndexRegion == r,1,'first');
            if ~isempty(index)
               rectangle(axRR,'Position',[0, (numRegions - r) * 1/numRegions, 1, 1/numRegions],'FaceColor',chanColors(index,:),'LineStyle','none')
            end
         end

         % Spectra Raw
         axSR = axes;
         axSR.Tag = 'axSR';
         axSR.Toolbar.Visible = 'off';
         axSR.Units = 'normalized';
         axSR.InnerPosition = [0,.45,1,.15];
         axSR.Box = 'off';
         axSR.Color = 'none';
         xLo = 0;
         xHi = max(freqs);
         yLo = nanmin(RTT_ColVec(spectraRaw(:,freqs <= floor(max(freqs) * .99))));
         yHi = nanmax(RTT_ColVec(spectraRaw(:,freqs >= 3)));
         [~,j] = sort(iqr(spectraRaw,2),'descend');
         hold on
         for k=1:size(spectraRaw,1)
            plot(axSR,freqs,spectraRaw(j(k),:),'Color',chanColors(j(k),:),'LineWidth',1)
         end
         plot(axSR,freqs,nanmedian(spectraRaw),'Color','y','LineWidth',4)
         for n=1:height(notchTable)
            if n == height(notchTable)
               line(axSR,[notchTable.Peak(n),notchTable.Peak(n)],[yLo,yHi * .05],'Color','r','LineWidth',1,'LineStyle','-')
               line(axSR,[notchTable.LoStop(n) + 1,notchTable.LoStop(n) + 1],[yLo,yHi * .05],'Color','k','LineWidth',1,'LineStyle',':')
               line(axSR,[notchTable.HiStop(n) - 1,notchTable.HiStop(n) - 1],[yLo,yHi * .05],'Color','k','LineWidth',1,'LineStyle',':')
            else
               line(axSR,[notchTable.Peak(n),notchTable.Peak(n)],[yLo,yHi],'Color','r','LineWidth',1,'LineStyle','-')
               line(axSR,[notchTable.LoStop(n) + 1,notchTable.LoStop(n) + 1],[yLo,yHi],'Color','k','LineWidth',1,'LineStyle',':')
               line(axSR,[notchTable.HiStop(n) - 1,notchTable.HiStop(n) - 1],[yLo,yHi],'Color','k','LineWidth',1,'LineStyle',':')
            end
         end
         axSR.XLim = [xLo,xHi];
         axSR.YLim = [yLo,yHi];
         axSR.XTick = [];
         axSR.YTick = [];
         axSR.XAxis.Visible = 'off';
         axSR.YAxis.Visible = 'off';

         % Spectra Notch
         axSN = axes;
         axSN.Tag = 'axSN';
         axSN.Toolbar.Visible = 'off';
         axSN.Units = 'normalized';
         axSN.InnerPosition = [0,.35,1,.15];
         axSN.Box = 'off';
         axSN.Color = 'none';
         xLo = 0;
         xHi = max(freqs);
         yLo = prctile(spectraNotch(:),1);
         yHi = nanmax(RTT_ColVec(spectraNotch(:,freqs >= 3)));
         [~,j] = sort(iqr(spectraNotch,2),'descend');
         hold on
         for j=1:size(spectraNotch,1)
            plot(axSN,freqs,spectraNotch(j,:),'Color',chanColors(j,:),'LineWidth',1)
         end   
         plot(axSN,freqs,nanmedian(spectraNotch),'Color','y','LineWidth',4)
         for n=1:height(notchTable)-1
            line(axSN,[notchTable.Peak(n),notchTable.Peak(n)],[yLo,yHi],'Color','r','LineWidth',1,'LineStyle','-')
            line(axSN,[notchTable.LoStop(n) + 1,notchTable.LoStop(n) + 1],[yLo,yHi],'Color','k','LineWidth',1,'LineStyle',':')
            line(axSN,[notchTable.HiStop(n) - 1,notchTable.HiStop(n) - 1],[yLo,yHi],'Color','k','LineWidth',1,'LineStyle',':')
         end
         axSN.XLim = [xLo,xHi];
         axSN.YLim = [yLo,yHi];
         axSN.XTick = [];
         axSN.YTick = [];
         axSN.XAxis.Visible = 'off';
         axSN.YAxis.Visible = 'off';

         % Spectra Physio
         axSP = axes;
         axSP.Tag = 'axSP';
         axSP.Toolbar.Visible = 'off';
         axSP.Units = 'normalized';
         axSP.InnerPosition = [0,.25,1,.15];
         axSP.Box = 'off';
         axSP.Color = 'none';
         xLo = 0;
         xHi = max(freqs);
         yLo = prctile(spectraPhysio(:),.1);
         yHi = nanmax(RTT_ColVec(spectraPhysio(:,freqs >= 3)));
         [~,j] = sort(iqr(spectraPhysio,2),'descend');
         hold on
         for j=1:size(spectraPhysio,1)
            plot(axSP,freqs,spectraPhysio(j,:),'Color',chanColors(j,:),'LineWidth',1)
         end   
         plot(axSP,freqs,nanmedian(spectraPhysio),'Color','y','LineWidth',4)
         for n=1:height(notchTable)
            line(axSP,[notchTable.Peak(n),notchTable.Peak(n)],[yLo,yHi],'Color','r','LineWidth',1,'LineStyle','-')
            line(axSP,[notchTable.LoStop(n) + 1,notchTable.LoStop(n) + 1],[yLo,yHi],'Color','k','LineWidth',1,'LineStyle',':')
            line(axSP,[notchTable.HiStop(n) - 1,notchTable.HiStop(n) - 1],[yLo,yHi],'Color','k','LineWidth',1,'LineStyle',':')
         end
         axSP.XLim = [xLo,xHi];
         axSP.YLim = [yLo,yHi];
         axSP.XTick = [];
         axSP.YTick = [];
         axSP.XAxis.Visible = 'off';
         axSP.YAxis.Visible = 'off';

         % Spectra Physio Inset
         axSPI = axes;
         axSPI.Tag = 'axSPI';
         axSPI.Toolbar.Visible = 'off';
         axSPI.Units = 'normalized';
         axSPI.InnerPosition = [0,0,.5,.25];
         axSPI.Box = 'off';
         axSPI.Color = 'none';
         spectraRegion = spectraPhysio(:,threshIndicesPhysio);
         xLo = 0;
         xHi = max(freqs(threshIndicesPhysio));
         yLo = nanmin(spectraRegion(:));
         yHi = nanmax(RTT_ColVec(spectraRegion(:,freqs(threshIndicesPhysio) >= 2)));
         [~,j] = sort(iqr(spectraRegion,2),'descend');
         hold on
         patch(axSPI,'Xdata',[threshBandDelta(1),threshBandDelta(2),threshBandDelta(2),threshBandDelta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','b','FaceAlpha',.15,'LineStyle','none')
         patch(axSPI,'Xdata',[threshBandTheta(1),threshBandTheta(2),threshBandTheta(2),threshBandTheta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','c','FaceAlpha',.15,'LineStyle','none')
         patch(axSPI,'Xdata',[threshBandAlpha(1),threshBandAlpha(2),threshBandAlpha(2),threshBandAlpha(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','g','FaceAlpha',.15,'LineStyle','none')
         patch(axSPI,'Xdata',[threshBandBeta1(1),threshBandBeta1(2),threshBandBeta1(2),threshBandBeta1(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','y','FaceAlpha',.15,'LineStyle','none')
         patch(axSPI,'Xdata',[threshBandBeta2(1),threshBandBeta2(2),threshBandBeta2(2),threshBandBeta2(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','m','FaceAlpha',.15,'LineStyle','none')
         patch(axSPI,'Xdata',[threshBandGamma(1),threshBandGamma(2),threshBandGamma(2),threshBandGamma(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','r','FaceAlpha',.15,'LineStyle','none')
         for j=1:size(spectraRegion,1)
            plot(axSPI,freqs(threshIndicesPhysio),spectraRegion(j,:),'Color',chanColors(j,:),'LineWidth',2)
         end   
         plot(axSPI,freqs(threshIndicesPhysio),nanmedian(spectraRegion),'Color','y','LineWidth',3)
         axSPI.XLim = [xLo,xHi];
         axSPI.YLim = [yLo,yHi];
         axSPI.XTick = [];
         axSPI.YTick = [];
         axSPI.XAxis.Visible = 'off';
         axSPI.YAxis.Visible = 'off';

         % Spectra Physio Region All
         axSPA = axes;
         axSPA.Tag = 'axSPA';
         axSPA.Toolbar.Visible = 'off';
         axSPA.Units = 'normalized';
         axSPA.InnerPosition = [.5,.125,.125,.125];
         axSPA.Box = 'off';
         axSPA.Color = 'none';
         spectraRegion = spectraPhysio(:,threshIndicesPhysio);
         xLo = 0;
         xHi = max(freqs(threshIndicesPhysio));
         yLo = nanmin(spectraRegion(:));
         yHi = nanmax(RTT_ColVec(spectraRegion(:,freqs(threshIndicesPhysio) >= 2)));
         [~,j] = sort(iqr(spectraRegion,2),'descend');
         hold on
         patch(axSPA,'Xdata',[threshBandDelta(1),threshBandDelta(2),threshBandDelta(2),threshBandDelta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','b','FaceAlpha',.15,'LineStyle','none')
         patch(axSPA,'Xdata',[threshBandTheta(1),threshBandTheta(2),threshBandTheta(2),threshBandTheta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','c','FaceAlpha',.15,'LineStyle','none')
         patch(axSPA,'Xdata',[threshBandAlpha(1),threshBandAlpha(2),threshBandAlpha(2),threshBandAlpha(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','g','FaceAlpha',.15,'LineStyle','none')
         patch(axSPA,'Xdata',[threshBandBeta1(1),threshBandBeta1(2),threshBandBeta1(2),threshBandBeta1(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','y','FaceAlpha',.15,'LineStyle','none')
         patch(axSPA,'Xdata',[threshBandBeta2(1),threshBandBeta2(2),threshBandBeta2(2),threshBandBeta2(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','m','FaceAlpha',.15,'LineStyle','none')
         patch(axSPA,'Xdata',[threshBandGamma(1),threshBandGamma(2),threshBandGamma(2),threshBandGamma(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','r','FaceAlpha',.15,'LineStyle','none')
         for j=1:size(spectraRegion,1)
            plot(axSPA,freqs(threshIndicesPhysio),spectraRegion(j,:),'Color','b','LineWidth',1)
         end   
         plot(axSPA,freqs(threshIndicesPhysio),nanmedian(spectraRegion),'Color','y','LineWidth',3)
         line(axSPA,[xLo,xHi],[yLo,yLo],'Color','k','LineWidth',2)
         axSPA.XLim = [xLo,xHi];
         axSPA.YLim = [yLo,yHi];
         axSPA.XTick = [];
         axSPA.YTick = [];
         axSPA.XAxis.Visible = 'off';
         axSPA.YAxis.Visible = 'off';

         % Spectra Physio Region Jaws
         r = 1;
         if ismember(r,chanlocsTable.IndexRegion)
            axSPJ = axes;
            axSPJ.Tag = 'axSPJ';
            axSPJ.UserData.IndexRegion = r;
            axSPJ.Toolbar.Visible = 'off';
            axSPJ.Units = 'normalized';
            axSPJ.InnerPosition = [.625,.125,.125,.125];
            axSPJ.Box = 'off';
            axSPJ.Color = 'none';
            spectraRegion = spectraPhysio(chanlocsTable.IndexRegion == r,threshIndicesPhysio);
            xLo = 0;
            xHi = max(freqs(threshIndicesPhysio));
            yLo = nanmin(spectraRegion(:));
            yHi = nanmax(RTT_ColVec(spectraRegion(:,freqs(threshIndicesPhysio) >= 2)));
            [~,j] = sort(iqr(spectraRegion,2),'descend');
            hold on
            patch(axSPJ,'Xdata',[threshBandDelta(1),threshBandDelta(2),threshBandDelta(2),threshBandDelta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','b','FaceAlpha',.15,'LineStyle','none')
            patch(axSPJ,'Xdata',[threshBandTheta(1),threshBandTheta(2),threshBandTheta(2),threshBandTheta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','c','FaceAlpha',.15,'LineStyle','none')
            patch(axSPJ,'Xdata',[threshBandAlpha(1),threshBandAlpha(2),threshBandAlpha(2),threshBandAlpha(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','g','FaceAlpha',.15,'LineStyle','none')
            patch(axSPJ,'Xdata',[threshBandBeta1(1),threshBandBeta1(2),threshBandBeta1(2),threshBandBeta1(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','y','FaceAlpha',.15,'LineStyle','none')
            patch(axSPJ,'Xdata',[threshBandBeta2(1),threshBandBeta2(2),threshBandBeta2(2),threshBandBeta2(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','m','FaceAlpha',.15,'LineStyle','none')
            patch(axSPJ,'Xdata',[threshBandGamma(1),threshBandGamma(2),threshBandGamma(2),threshBandGamma(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','r','FaceAlpha',.15,'LineStyle','none')
            for j=1:size(spectraRegion,1)
               plot(axSPJ,freqs(threshIndicesPhysio),spectraRegion(j,:),'LineWidth',1)
            end   
            plot(axSPJ,freqs(threshIndicesPhysio),nanmedian(spectraRegion),'Color',chanColorsRegions(r,:),'LineWidth',3)
            line(axSPJ,[xLo,xHi],[yLo,yLo],'Color','k','LineWidth',2)
            if isnan(yLo)
               yLo = 0;
               yHi = 1;
            else
               axSPJ.XLim = [xLo,xHi];
               axSPJ.YLim = [yLo,yHi];
            end
            axSPJ.XTick = [];
            axSPJ.YTick = [];
            axSPJ.XAxis.Visible = 'off';
            axSPJ.YAxis.Visible = 'off';
         end

         % Spectra Physio Region Eyes
         r = 2;
         axSPE = axes;
         axSPE.Tag = 'axSPE';
         axSPE.UserData.IndexRegion = r;
         axSPE.Toolbar.Visible = 'off';
         axSPE.Units = 'normalized';
         axSPE.InnerPosition = [.75,.125,.125,.125];
         axSPE.Box = 'off';
         axSPE.Color = 'none';
         spectraRegion = spectraPhysio(chanlocsTable.IndexRegion == r,threshIndicesPhysio);
         xLo = 0;
         xHi = max(freqs(threshIndicesPhysio));
         yLo = nanmin(spectraRegion(:));
         yHi = nanmax(RTT_ColVec(spectraRegion(:,freqs(threshIndicesPhysio) >= 2)));
         [~,j] = sort(iqr(spectraRegion,2),'descend');
         hold on
         patch(axSPE,'Xdata',[threshBandDelta(1),threshBandDelta(2),threshBandDelta(2),threshBandDelta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','b','FaceAlpha',.15,'LineStyle','none')
         patch(axSPE,'Xdata',[threshBandTheta(1),threshBandTheta(2),threshBandTheta(2),threshBandTheta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','c','FaceAlpha',.15,'LineStyle','none')
         patch(axSPE,'Xdata',[threshBandAlpha(1),threshBandAlpha(2),threshBandAlpha(2),threshBandAlpha(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','g','FaceAlpha',.15,'LineStyle','none')
         patch(axSPE,'Xdata',[threshBandBeta1(1),threshBandBeta1(2),threshBandBeta1(2),threshBandBeta1(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','y','FaceAlpha',.15,'LineStyle','none')
         patch(axSPE,'Xdata',[threshBandBeta2(1),threshBandBeta2(2),threshBandBeta2(2),threshBandBeta2(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','m','FaceAlpha',.15,'LineStyle','none')
         patch(axSPE,'Xdata',[threshBandGamma(1),threshBandGamma(2),threshBandGamma(2),threshBandGamma(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','r','FaceAlpha',.15,'LineStyle','none')
         for j=1:size(spectraRegion,1)
            plot(axSPE,freqs(threshIndicesPhysio),spectraRegion(j,:),'LineWidth',1)
         end   
         plot(axSPE,freqs(threshIndicesPhysio),nanmedian(spectraRegion),'Color',chanColorsRegions(r,:),'LineWidth',3)
         line(axSPE,[xLo,xHi],[yLo,yLo],'Color','k','LineWidth',2)
         if isnan(yLo)
            yLo = 0;
            yHi = 1;
         else
            axSPE.XLim = [xLo,xHi];
            axSPE.YLim = [yLo,yHi];
         end
         axSPE.XTick = [];
         axSPE.YTick = [];
         axSPE.XAxis.Visible = 'off';
         axSPE.YAxis.Visible = 'off';

         % Spectra Physio Region Forehead
         r = 3;
         axSPF = axes;
         axSPF.Tag = 'axSPF';
         axSPF.UserData.IndexRegion = r;
         axSPF.Toolbar.Visible = 'off';
         axSPF.Units = 'normalized';
         axSPF.InnerPosition = [.875,.125,.125,.125];
         axSPF.Box = 'off';
         axSPF.Color = 'none';
         spectraRegion = spectraPhysio(chanlocsTable.IndexRegion == r,threshIndicesPhysio);
         xLo = 0;
         xHi = max(freqs(threshIndicesPhysio));
         yLo = nanmin(spectraRegion(:));
         yHi = nanmax(RTT_ColVec(spectraRegion(:,freqs(threshIndicesPhysio) >= 2)));
         [~,j] = sort(iqr(spectraRegion,2),'descend');
         hold on
         patch(axSPF,'Xdata',[threshBandDelta(1),threshBandDelta(2),threshBandDelta(2),threshBandDelta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','b','FaceAlpha',.15,'LineStyle','none')
         patch(axSPF,'Xdata',[threshBandTheta(1),threshBandTheta(2),threshBandTheta(2),threshBandTheta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','c','FaceAlpha',.15,'LineStyle','none')
         patch(axSPF,'Xdata',[threshBandAlpha(1),threshBandAlpha(2),threshBandAlpha(2),threshBandAlpha(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','g','FaceAlpha',.15,'LineStyle','none')
         patch(axSPF,'Xdata',[threshBandBeta1(1),threshBandBeta1(2),threshBandBeta1(2),threshBandBeta1(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','y','FaceAlpha',.15,'LineStyle','none')
         patch(axSPF,'Xdata',[threshBandBeta2(1),threshBandBeta2(2),threshBandBeta2(2),threshBandBeta2(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','m','FaceAlpha',.15,'LineStyle','none')
         patch(axSPF,'Xdata',[threshBandGamma(1),threshBandGamma(2),threshBandGamma(2),threshBandGamma(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','r','FaceAlpha',.15,'LineStyle','none')
         for j=1:size(spectraRegion,1)
            plot(axSPF,freqs(threshIndicesPhysio),spectraRegion(j,:),'LineWidth',1)
         end   
         plot(axSPF,freqs(threshIndicesPhysio),nanmedian(spectraRegion),'Color',chanColorsRegions(r,:),'LineWidth',3)
         line(axSPF,[xLo,xHi],[yLo,yLo],'Color','k','LineWidth',2)
         if isnan(yLo)
            yLo = 0;
            yHi = 1;
         else
            axSPF.XLim = [xLo,xHi];
            axSPF.YLim = [yLo,yHi];
         end
         axSPF.XTick = [];
         axSPF.YTick = [];
         axSPF.XAxis.Visible = 'off';
         axSPF.YAxis.Visible = 'off';

         % Spectra Physio Region Top
         r = 4;
         axSPT = axes;
         axSPT.Tag = 'axSPT';
         axSPT.UserData.IndexRegion = r;
         axSPT.Toolbar.Visible = 'off';
         axSPT.Units = 'normalized';
         axSPT.InnerPosition = [.5,0,.125,.125];
         axSPT.Box = 'off';
         axSPT.Color = 'none';
         spectraRegion = spectraPhysio(chanlocsTable.IndexRegion == r,threshIndicesPhysio);
         xLo = 0;
         xHi = max(freqs(threshIndicesPhysio));
         yLo = nanmin(spectraRegion(:));
         yHi = nanmax(RTT_ColVec(spectraRegion(:,freqs(threshIndicesPhysio) >= 2)));
         [~,j] = sort(iqr(spectraRegion,2),'descend');
         hold on
         patch(axSPT,'Xdata',[threshBandDelta(1),threshBandDelta(2),threshBandDelta(2),threshBandDelta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','b','FaceAlpha',.15,'LineStyle','none')
         patch(axSPT,'Xdata',[threshBandTheta(1),threshBandTheta(2),threshBandTheta(2),threshBandTheta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','c','FaceAlpha',.15,'LineStyle','none')
         patch(axSPT,'Xdata',[threshBandAlpha(1),threshBandAlpha(2),threshBandAlpha(2),threshBandAlpha(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','g','FaceAlpha',.15,'LineStyle','none')
         patch(axSPT,'Xdata',[threshBandBeta1(1),threshBandBeta1(2),threshBandBeta1(2),threshBandBeta1(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','y','FaceAlpha',.15,'LineStyle','none')
         patch(axSPT,'Xdata',[threshBandBeta2(1),threshBandBeta2(2),threshBandBeta2(2),threshBandBeta2(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','m','FaceAlpha',.15,'LineStyle','none')
         patch(axSPT,'Xdata',[threshBandGamma(1),threshBandGamma(2),threshBandGamma(2),threshBandGamma(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','r','FaceAlpha',.15,'LineStyle','none')
         for j=1:size(spectraRegion,1)
            plot(axSPT,freqs(threshIndicesPhysio),spectraRegion(j,:),'LineWidth',1)
         end   
         plot(axSPT,freqs(threshIndicesPhysio),nanmedian(spectraRegion),'Color',chanColorsRegions(r,:),'LineWidth',3)
         if isnan(yLo)
            yLo = 0;
            yHi = 1;
         else
            axSPT.XLim = [xLo,xHi];
            axSPT.YLim = [yLo,yHi];
         end
         axSPT.XTick = [];
         axSPT.YTick = [];
         axSPT.XAxis.Visible = 'off';
         axSPT.YAxis.Visible = 'off';

         % Spectra Physio Region Left
         r = 5;
         axSPL = axes;
         axSPL.Tag = 'axSPL';
         axSPL.UserData.IndexRegion = r;
         axSPL.Toolbar.Visible = 'off';
         axSPL.Units = 'normalized';
         axSPL.InnerPosition = [.625,0,.125,.125];
         axSPL.Box = 'off';
         axSPL.Color = 'none';
         spectraRegion = spectraPhysio(chanlocsTable.IndexRegion == r,threshIndicesPhysio);
         xLo = 0;
         xHi = max(freqs(threshIndicesPhysio));
         yLo = nanmin(spectraRegion(:));
         yHi = nanmax(RTT_ColVec(spectraRegion(:,freqs(threshIndicesPhysio) >= 2)));
         [~,j] = sort(iqr(spectraRegion,2),'descend');
         hold on
         patch(axSPL,'Xdata',[threshBandDelta(1),threshBandDelta(2),threshBandDelta(2),threshBandDelta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','b','FaceAlpha',.15,'LineStyle','none')
         patch(axSPL,'Xdata',[threshBandTheta(1),threshBandTheta(2),threshBandTheta(2),threshBandTheta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','c','FaceAlpha',.15,'LineStyle','none')
         patch(axSPL,'Xdata',[threshBandAlpha(1),threshBandAlpha(2),threshBandAlpha(2),threshBandAlpha(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','g','FaceAlpha',.15,'LineStyle','none')
         patch(axSPL,'Xdata',[threshBandBeta1(1),threshBandBeta1(2),threshBandBeta1(2),threshBandBeta1(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','y','FaceAlpha',.15,'LineStyle','none')
         patch(axSPL,'Xdata',[threshBandBeta2(1),threshBandBeta2(2),threshBandBeta2(2),threshBandBeta2(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','m','FaceAlpha',.15,'LineStyle','none')
         patch(axSPL,'Xdata',[threshBandGamma(1),threshBandGamma(2),threshBandGamma(2),threshBandGamma(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','r','FaceAlpha',.15,'LineStyle','none')
         for j=1:size(spectraRegion,1)
            plot(axSPL,freqs(threshIndicesPhysio),spectraRegion(j,:),'LineWidth',1)
         end   
         plot(axSPL,freqs(threshIndicesPhysio),nanmedian(spectraRegion),'Color',chanColorsRegions(r,:),'LineWidth',3)
         if isnan(yLo)
            yLo = 0;
            yHi = 1;
         else
            axSPL.XLim = [xLo,xHi];
            axSPL.YLim = [yLo,yHi];
         end
         axSPL.XTick = [];
         axSPL.YTick = [];
         axSPL.XAxis.Visible = 'off';
         axSPL.YAxis.Visible = 'off';

         % Spectra Physio Region Right
         r = 6;
         axSPR = axes;
         axSPR.Tag = 'axSPR';
         axSPR.UserData.IndexRegion = r;
         axSPR.Toolbar.Visible = 'off';
         axSPR.Units = 'normalized';
         axSPR.InnerPosition = [.75,0,.125,.125];
         axSPR.Box = 'off';
         axSPR.Color = 'none';
         spectraRegion = spectraPhysio(chanlocsTable.IndexRegion == r,threshIndicesPhysio);
         xLo = 0;
         xHi = max(freqs(threshIndicesPhysio));
         yLo = nanmin(spectraRegion(:));
         yHi = nanmax(RTT_ColVec(spectraRegion(:,freqs(threshIndicesPhysio) >= 2)));
         [~,j] = sort(iqr(spectraRegion,2),'descend');
         hold on
         patch(axSPR,'Xdata',[threshBandDelta(1),threshBandDelta(2),threshBandDelta(2),threshBandDelta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','b','FaceAlpha',.15,'LineStyle','none')
         patch(axSPR,'Xdata',[threshBandTheta(1),threshBandTheta(2),threshBandTheta(2),threshBandTheta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','c','FaceAlpha',.15,'LineStyle','none')
         patch(axSPR,'Xdata',[threshBandAlpha(1),threshBandAlpha(2),threshBandAlpha(2),threshBandAlpha(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','g','FaceAlpha',.15,'LineStyle','none')
         patch(axSPR,'Xdata',[threshBandBeta1(1),threshBandBeta1(2),threshBandBeta1(2),threshBandBeta1(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','y','FaceAlpha',.15,'LineStyle','none')
         patch(axSPR,'Xdata',[threshBandBeta2(1),threshBandBeta2(2),threshBandBeta2(2),threshBandBeta2(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','m','FaceAlpha',.15,'LineStyle','none')
         patch(axSPR,'Xdata',[threshBandGamma(1),threshBandGamma(2),threshBandGamma(2),threshBandGamma(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','r','FaceAlpha',.15,'LineStyle','none')
         for j=1:size(spectraRegion,1)
            plot(axSPR,freqs(threshIndicesPhysio),spectraRegion(j,:),'LineWidth',1)
         end   
         plot(axSPR,freqs(threshIndicesPhysio),nanmedian(spectraRegion),'Color',chanColorsRegions(r,:),'LineWidth',3)
         if isnan(yLo)
            yLo = 0;
            yHi = 1;
         else
            axSPR.XLim = [xLo,xHi];
            axSPR.YLim = [yLo,yHi];
         end
         axSPR.XTick = [];
         axSPR.YTick = [];
         axSPR.XAxis.Visible = 'off';
         axSPR.YAxis.Visible = 'off';

         % Spectra Physio Region Back
         r = numRegions;
         axSPB = axes;
         axSPB.Tag = 'axSPB';
         axSPB.UserData.IndexRegion = r;
         axSPB.Toolbar.Visible = 'off';
         axSPB.Units = 'normalized';
         axSPB.InnerPosition = [.875,0,.125,.125];
         axSPB.Box = 'off';
         axSPB.Color = 'none';
         spectraRegion = spectraPhysio(chanlocsTable.IndexRegion == r,threshIndicesPhysio);
         xLo = 0;
         xHi = max(freqs(threshIndicesPhysio));
         yLo = nanmin(spectraRegion(:));
         yHi = nanmax(RTT_ColVec(spectraRegion(:,freqs(threshIndicesPhysio) >= 2)));
         [~,j] = sort(iqr(spectraRegion,2),'descend');
         hold on
         patch(axSPB,'Xdata',[threshBandDelta(1),threshBandDelta(2),threshBandDelta(2),threshBandDelta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','b','FaceAlpha',.15,'LineStyle','none')
         patch(axSPB,'Xdata',[threshBandTheta(1),threshBandTheta(2),threshBandTheta(2),threshBandTheta(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','c','FaceAlpha',.15,'LineStyle','none')
         patch(axSPB,'Xdata',[threshBandAlpha(1),threshBandAlpha(2),threshBandAlpha(2),threshBandAlpha(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','g','FaceAlpha',.15,'LineStyle','none')
         patch(axSPB,'Xdata',[threshBandBeta1(1),threshBandBeta1(2),threshBandBeta1(2),threshBandBeta1(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','y','FaceAlpha',.15,'LineStyle','none')
         patch(axSPB,'Xdata',[threshBandBeta2(1),threshBandBeta2(2),threshBandBeta2(2),threshBandBeta2(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','m','FaceAlpha',.15,'LineStyle','none')
         patch(axSPB,'Xdata',[threshBandGamma(1),threshBandGamma(2),threshBandGamma(2),threshBandGamma(1)],'YData',[yLo,yLo,yHi,yHi],'FaceColor','r','FaceAlpha',.15,'LineStyle','none')
         for j=1:size(spectraRegion,1)
            plot(axSPB,freqs(threshIndicesPhysio),spectraRegion(j,:),'LineWidth',1)
         end   
         plot(axSPB,freqs(threshIndicesPhysio),nanmedian(spectraRegion),'Color',chanColorsRegions(r,:),'LineWidth',3)
         if isnan(yLo)
            yLo = 0;
            yHi = 1;
         else
            axSPB.XLim = [xLo,xHi];
            axSPB.YLim = [yLo,yHi];
         end
         axSPB.XTick = [];
         axSPB.YTick = [];
         axSPB.XAxis.Visible = 'off';
         axSPB.YAxis.Visible = 'off';

         % Interface
         chansBadIter = [];
         dataRBprev = dataRB;
         spectraRawPrev = spectraRaw;
         spectraNotchPrev = spectraNotch;
         spectraPhysioPrev = spectraPhysio;

         repeat = true;
         while repeat
            txSPVP.Color = 'r';
            indexFrequency = [];
            c = [];
            [mouseX,mouseY,mouseB] = ginput(1);
            axCurr = fig.CurrentAxes;
            if mouseB == 3
               % undo
               txSPVP.Color = 'y';
               pause(.1)
               repeat = false;
               recompute = true;
               chansBadIter = [];
               dataRB = dataRBprev;
               spectraRaw = spectraRawPrev;
               spectraNotch = spectraNotchPrev;
               spectraPhysio = spectraPhysioPrev;
            elseif mouseB == 1
               if startsWith(axCurr.Tag,'axS')
                  % get frequency index
                  [~,indexFrequency] = min(abs(mouseX - freqs));
                  if isfield(axCurr.UserData,'IndexRegion')
                     spectraRegion = spectraPhysio;
                     spectraRegion(chanlocsTable.IndexRegion ~= axCurr.UserData.IndexRegion,:) = NaN;
                     if all(isnan(RTT_ColVec(spectraRegion(chanlocsTable.IndexRegion == axCurr.UserData.IndexRegion,:))))
                        c = find(chanlocsTable.IndexRegion == axCurr.UserData.IndexRegion,1,'first');
                     else
                        [~,c] = nanmin(abs(mouseY - spectraRegion(:,indexFrequency)));
                     end
                  elseif strcmp(axCurr.Tag,'axSR')
                     [~,c] = nanmin(abs(mouseY - spectraRaw(:,indexFrequency)));
                  elseif strcmp(axCurr.Tag,'axSN')
                     [~,c] = nanmin(abs(mouseY - spectraNotch(:,indexFrequency)));
                  else
                     [~,c] = nanmin(abs(mouseY - spectraPhysio(:,indexFrequency)));         
                  end
                  % update chandBadIter
                  chansBadIter = unique([chansBadIter ; c]);
                  % white out channels
                  plot(axRB,t,dataRB(c,:),'w','LineWidth',1)
                  plot(axSR,freqs,spectraRaw(c,:),'w','LineWidth',1)
                  plot(axSN,freqs,spectraNotch(c,:),'w','LineWidth',1)
                  plot(axSP,freqs,spectraPhysio(c,:),'w','LineWidth',1)
                  plot(axSPI,freqs,spectraPhysio(c,:),'w','LineWidth',2)
                  plot(axSPA,freqs,spectraPhysio(c,:),'w','LineWidth',1)
                  eval(sprintf('axPlot = axSP%c;',chanlocsTable.Region{c}(1)))
                  plot(axPlot,freqs,spectraPhysio(c,:),'w','LineWidth',1)
                  % update dataRB and spectra
                  dataRB(chansBadIter,:) = NaN;
                  spectraRaw(chansBadIter,:) = NaN;
                  spectraNotch(chansBadIter,:) = NaN;
                  spectraPhysio(chansBadIter,:) = NaN;
               elseif isempty(chansBadIter)
                  % bad channels complete, manual renotching decision
                  repeat = false;
                  txSPVP.String = strcat({'RIGHT CLICK TO FLAG???   '},threshFigTitle);
                  txSPVP.Color = 'm';
                  pause(.25)
                  if isempty(chansBad)
                     txBC.String = sprintf('Bad Chans: %u (%g%%)\n',length(chansBad),round(100 * length(chansBad) / nbchan,1));
                  else
                     txBC.String = sprintf('Bad Chans: %u (%g%%)\n%s',length(chansBad),round(100 * length(chansBad) / nbchan,1),strtrim(regexprep(mat2str(chansBad),'\D','\n')));
                  end
                  clearvars *Prev spectraRegion
                  [~,~,mouseB] = ginput(1);
                  if mouseB == 1
                     % complete
                     txSPVP.String = threshFigTitle;
                     txSPVP.Color = 'k';
                     pause(.25)
                     % save
                     chanlocsTable.IsBad(chansBad) = true;
                     chanlocsTable.ReasonBad = repmat({''},[height(chanlocsTable),1]);
                     chanlocsTable.ReasonBad(chansBad) = {'spectra'};
                  else
                     % flag output
                     txSPVP.String = strcat({'   !!! INSPECT FLAG !!!   '},threshFigTitle,{'   !!! INSPECT FLAG !!!   '});
                     txSPVP.Color = 'r';
                     D.PathAndNameNewInspected{iFile} = strrep(D.PathAndNameNewInspected{iFile},'_INSPECTED.mat','_!REINSPECT.mat');
                     D.PathAndNameNewFigure{iFile} = strrep(D.PathAndNameNewFigure{iFile},'_INSPECTED.jpg','_!REINSPECT.jpg');
                     pause(.25)
                  end
                  save(D.PathAndNameNewInspected{iFile},'spectra*','freqs','*Table','t','data')
               else
                  txSPVP.Color = 'y';
                  pause(.1)
                  repeat = false;
                  recompute = true;
                  chansBad = unique([chansBad ; chansBadIter]);
               end
            else
               error('unrecognized mouse button')
            end
            % status update
            fprintf('%s, f = %u Hz, dB = %g, c = %u, badIter = %s, chansBad = %s\n',axCurr.Tag,freqs(indexFrequency),round(mouseY,2),c,mat2str(chansBadIter),mat2str(chansBad))
         end
      end

      % Save Figure
      exportgraphics(fig,D.PathAndNameNewFigure{iFile},'Resolution',threshFigRes)
      close(fig)
   end
end
disp('section 4 complete')

%% 5. Filter, interpolate bad channels, rereference, resample, run ICA

clearvars
clc

folderPreProc2 = 'R:\COLLAB\BABY\EEG\PREPROC\PREPROC2';
D = RTT_DirList(folderPreProc2);
D = D(strcmp(D.Ext,'.set'),:);
D.PathAndNameNew = strrep(D.PathAndName,'PROC2','PROC3');
D.PathAndNameNewInspected = strrep(D.PathAndName,'.set','_INSPECTED.mat');

for iFile=1:height(D)
   if ~exist(D.PathAndNameNew{iFile},'file') && exist(D.PathAndNameNewInspected{iFile},'file')
      EEG = pop_loadset(D.PathAndName{iFile});
      load(D.PathAndNameNewInspected{iFile})
      % filter
      EEG = pop_eegfiltnew(EEG,2,[]); % cutoff 1 hz
      EEG = pop_eegfiltnew(EEG,[],40); % cutoff 45 hz
      EEG.chanlocs = table2struct(chanlocsTable);
      if sum(chanlocsTable.IsBad) > 0
         EEG = pop_interp(EEG,find(chanlocsTable.IsBad),'spherical');
      end
      % limit data
      % extreme magnitude is limited to improve the ICA algorithm's efficiency
      threshMag = 200; % uV
      EEG.data(EEG.data > threshMag) = threshMag;
      EEG.data(EEG.data < -threshMag) = -threshMag;
      % reref
      EEG = pop_reref(EEG,[]);
      % resamp
      if EEG.srate > 256
         EEG = pop_resample(EEG,250);
      end
      % run infomax ica
      EEG = pop_runica(EEG,'extended',1,'pca',rank(EEG.data) - 2);
      % generate automatic component labels
      EEG = iclabel(EEG,'default');
      % save
      EEG = pop_saveset(EEG,'filename',D.PathAndNameNew{iFile},'savemode','onefile','version','7.3');
   end
end
disp('section 5 complete')

%% 6. Manually inspect ICA decompositions and save .mat with bad components identified for pruning

% ICA decomposition is inspected manually by a human expert, bad components are identified and saved as a .mat in the
% same folder

clearvars
clc

folderProc3 = 'R:\COLLAB\BABY\EEG\PREPROC\PREPROC3';
D = RTT_DirList(folderProc3);
D = D(strcmp(D.Ext,'.set'),:);
D.PathAndNameNewPruned = strrep(D.PathAndName,'.set','_PRUNED.mat');

for iFile=1:height(D)
   if ~exist(D.PathAndNameNewPruned{iFile},'file')
      load('R:\COLLAB\BABY\EEG\ADMIN\masterEmptyGlobalsEEGLAB.mat')
      EEG = pop_loadset(D.PathAndName{iFile});
      EEG = pop_selectcomps(EEG,1:15);
      fig = gcf;
      fig.WindowState = 'maximized';
      waitfor(fig)
      compsPruned = EEG.reject.gcompreject;
      save(D.PathAndNameNewPruned{iFile},'compsPruned');
      clearvars -except D iFile
   end
end

disp('section 6 complete')

%% 7. Subtract artifactual components, rereference and save the final clean .set

clearvars
clc

folderProc3 = 'R:\COLLAB\BABY\EEG\PREPROC\PREPROC3';
D = RTT_DirList(folderProc3);
D = D(strcmp(D.Ext,'.set'),:);
D.PathAndNameNew = strrep(D.PathAndName,'PROC3','FINAL');
D.PathAndNameNewPruned = strrep(D.PathAndName,'.set','_PRUNED.mat');

for iFile=1:height(D)
   if exist(D.PathAndNameNewPruned{iFile},'file') && ~exist(D.PathAndNameNew{iFile},'file')
      EEG = pop_loadset(D.PathAndName{iFile});
      load(D.PathAndNameNewPruned{iFile})
      if sum(compsPruned) > 0
         EEG = pop_subcomp(EEG,find(compsPruned));
         % reref
         EEG = pop_reref(EEG,[]);         
      end
      EEG = pop_saveset(EEG,'filename',D.PathAndNameNew{iFile},'savemode','onefile','version','7.3');
   end
end

disp('section 7 complete')

%% custom functions

function columnVector = RTT_ColVec(inputArrayOrCell)
   if isnumeric(inputArrayOrCell) || islogical(inputArrayOrCell)
      columnVector = inputArrayOrCell(:);
   elseif iscell(inputArrayOrCell)
      columnVector = vertcat(inputArrayOrCell{:});
   else
      error('invalid data type')
   end
end

function D = RTT_DirList(folderOrFiles)
   if ~iscell(folderOrFiles)
      [~,folderOrFiles] = RTT_FoldersAndFilesRecursive(folderOrFiles);
   end
   D = cellfun(@(x) (dir(x)),folderOrFiles,'UniformOutput',false);
   D = vertcat(D{:});
   D = struct2table(D);
   D.name = cellstr(D.name);
   D.folder = cellstr(D.folder);
   D.date = cellstr(D.date);
   D = table(D.folder,D.name,D.bytes,D.datenum,'VariableNames',{'Path','Name','Bytes','DN'});
   D.PathAndName = fullfile(D.Path,D.Name);
   [~,~,D.Ext] = cellfun(@fileparts,D.PathAndName,'UniformOutput',false);
   D.Name = cellfun(@(x,y) x(1:end-length(y)),D.Name,D.Ext,'UniformOutput',false); % updated
   D.TS = RTT_Datestr(D.DN);
   D.Index = (1:height(D))';
   D = D(:,{'Index','Path','Name','Ext','PathAndName','Bytes','DN','TS'});
end

function [allFolders,allFullPaths] = RTT_FoldersAndFilesRecursive(startFolder)
   allFolders = cell('');
   allFullPaths = cell('');
   currentFolders = {startFolder};
   keepGoing = true;
   while keepGoing
      if ~isempty(currentFolders)
         allFolders = vertcat(allFolders,currentFolders);
         dirOutput = cellfun(@dir,currentFolders,'UniformOutput',false);

         filesOutput = cellfun(@(x) x(~[x.isdir]),dirOutput,'UniformOutput',false);
         filesOutputFullPath = cellfun(@(x) (fullfile({x.folder}',{x.name}')),filesOutput,'UniformOutput',false);
         filesOutputFullPathConcat = vertcat(filesOutputFullPath{:});
         allFullPaths = vertcat(allFullPaths,filesOutputFullPathConcat);

         folderOutput = cellfun(@(x) x(~ismember({x.name},{'.','..'}) & [x.isdir]),dirOutput,'UniformOutput',false);
         folderOutputFullPath = cellfun(@(x) (fullfile({x.folder}',{x.name}')),folderOutput,'UniformOutput',false);
         folderOutputFullPathConcat = vertcat(folderOutputFullPath{:});
         if isempty(folderOutputFullPathConcat)
            keepGoing = false;
         else
            currentFolders = folderOutputFullPathConcat;
         end
      end
   end
   allFolders = allFolders(~strcmp(allFolders,startFolder));
end

function ts = RTT_Datestr(dn)
   if isnumeric(dn)
      dn(isnan(dn)) = 0;
      ts = upper(datestr(dn,'yyyy-mm-dd_HHMM_ddd'));
   elseif isa(dn,'cell')
      ts = upper(cellfun(@(x) datestr(x,'yyyy-mm-dd_HHMM_ddd'),dn,'UniformOutput',false));
   end
   if ischar(ts)
      ts = cellstr(ts);
   end
   ts(startsWith(ts,'0000')) = {''};
   if all(contains(ts,'_0000_'))
      ts = erase(ts,'_0000');
   end
end
