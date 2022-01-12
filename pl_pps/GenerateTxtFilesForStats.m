function [fcmatrix_all] = GenerateTxtFilesForStats(cfg);
% This program will combine the fc values between/within the major lobes, separately for different frequency bands;
if ~isfield(cfg, 'foi'),              cfg.foi              = 'theta';       end; foi        = cfg.foi;
if ~isfield(cfg, 'fcmethod'),         cfg.fcmethod         = 'AAC';         end; fcmethod   = cfg.fcmethod;
if ~isfield(cfg, 'replacearg'),       cfg.replacearg       = 0;             end; replacearg = cfg.replacearg;
if ~isfield(cfg, 'parcmethod'),       cfg.parcmethod       = 'average';     end; parcmethod = cfg.parcmethod;
if ~isfield(cfg, 'methodtype'),       cfg.methodtype       = 'eloreta';     end; methodtype = cfg.methodtype;
if ~isfield(cfg, 'atlastype'),        cfg.atlastype        = 'LPBA';        end; atlastype  = cfg.atlastype;
if ~isfield(cfg, 'gridresolution'),   cfg.gridresolution   = '6mm';         end; gridresolution     = cfg.gridresolution;
if ~isfield(cfg, 'age'),              cfg.age              = 12;            end; age        = cfg.age;
%Check file if exists;
ouputfolder ='../Outputs/';
%outputname = ['Participant ' num2str(participantnumber) ' Age ' num2str(age) ' fcanalysis_' methodtype '_' foi '_' parcmethod '_' gridresolution '_' fcmethod '_' atlastype '.mat'];

global fcmatrix fcvalues_all;

%% get the participantnumbers
%ages = [12 36];
participantnumbers  = [];
%for i = 1:length(ages)
%age = ages(i);
datapath = ['../Data/RAW/Age' num2str(age) 'mos/'];
dirinfo = dir(datapath);
%dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
ii = 0;
for K = 1 : length(dirinfo)
if strcmp(dirinfo(K).name,'.') | strcmp(dirinfo(K).name,'..')
else
    ii = ii + 1;
    thisdir = dirinfo(K).name;
    thisstring = extractBetween(thisdir,"Subject "," Age");
        subjects(ii) = str2num(thisstring{1});
    end
end
participantnumbers.(['age' num2str(age)]) = unique(subjects);
%end 


nsubjs = length(participantnumbers.(['age' num2str(age)]));
disp(['The following analysis will be run for ' num2str(nsubjs) ' Participants.']);


fcmatrix_all  = [];
fcvalues_all  = zeros(nsubjs,11);

%% Get the data from each participant;
for jj = 1:nsubjs
participantnumber = participantnumbers.(['age' num2str(age)])(jj);
disp(['Participantnumber ' num2str(participantnumber) ' Age ' num2str(age)]);
% define the frequency band
switch foi % change the bpfreq boundaries according to age
    case 'theta'
        if age == 12
            bpfreq = [3 6]; 
        else % 36mos
            bpfreq = [3 7]; 
        end      
    case 'alpha'
        if age == 12
            bpfreq = [5 10];
        else
            bpfreq = [6 11];
        end
    case 'beta'
            bpfreq = [11 22];
    case 'gamma'
            bpfreq = [22 45];
end

% get the overall fc connectivity value from the LPBA atlas (48 Cortical ROIs);
datapth    = '../Outputs/';
matrixname = ['Participant ' num2str(participantnumber) ' Age ' num2str(age) ' fcanalysis_' methodtype '_' foi '_' parcmethod '_' gridresolution '_' fcmethod '_' atlastype '.mat']; 
load([datapth matrixname], 'fcmatrix');
nroi = 48; %48 cortical regions in the LPBA atlas;
switch fcmethod
    case 'AAC'
        junk = fcmatrix.Ro(1:nroi,1:nroi);
    case 'wpli'
        junk = fcmatrix.wpli(1:nroi,1:nroi);
    case {'coh','imag'}
        junk = fcmatrix.coh(1:nroi,1:nroi);
end
junk(isnan(junk))          = 0;
junk  = zscore(junk);
thres = prctile(junk(:),50); % top 50% of the data;
junk(junk<thres) = 0;
numofconns  = (nroi^2-nroi)*0.5;  %50% of the values were added into the sum value;
lpbafcvalues = sum(junk(:))/numofconns;

fcmatrix_all(jj,:,:) = junk;
%% Get the functional connectivity values between/within the four major lobes. 
 % connections within and between the major lobes;
lobarconn = zeros(1,10);
labelnums = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;3;3;3;3;3;3;3;3;3;3;4;4;4;4;4;4;4;4;2;2;2;2;2;2;2;2;4;4;2;2;2;2;3;3]; % WX made these indexes from the LPBA40 atlas labels;
% size_m = size(junk);
fidx   = find(labelnums ==1)';   
tidx   = find(labelnums ==2)';
pidx   = find(labelnums ==3)';
oidx   = find(labelnums ==4)';
%
switch fcmethod
    case 'AAC'
        junk = fcmatrix.Ro(1:nroi,1:nroi);
    case 'wpli'
        junk = fcmatrix.wpli(1:nroi,1:nroi);
    case {'coh','imag'}
        junk = fcmatrix.coh(1:nroi,1:nroi);
end
junk(isnan(junk))          = 0;
thres = prctile(junk(:),50);  
junk(junk<thres) = 0;
junk = zscore(junk); % take the z-score;
%
junk1 = junk(fidx,fidx);
lobarconn(1) = sum(junk1(:)) ./ length(find(junk1)); % Take the average of the nonzero connections;
junk1 = junk(fidx,tidx);
lobarconn(2) = sum(junk1(:)) ./ length(find(junk1));
junk1 = junk(fidx,pidx);
lobarconn(3) = sum(junk1(:)) ./ length(find(junk1));
junk1 = junk(fidx,oidx);
lobarconn(4) = sum(junk1(:)) ./ length(find(junk1));
junk1 = junk(tidx,tidx);
lobarconn(5) = sum(junk1(:)) ./ length(find(junk1));
junk1 = junk(tidx,pidx);
lobarconn(6) = sum(junk1(:)) ./ length(find(junk1));
junk1 = junk(tidx,oidx);
lobarconn(7) = sum(junk1(:)) ./ length(find(junk1));
junk1 = junk(pidx,pidx);
lobarconn(8) = sum(junk1(:)) ./ length(find(junk1));
junk1 = junk(pidx,oidx);
lobarconn(9) = sum(junk1(:)) ./ length(find(junk1));
junk1 = junk(oidx,oidx);
lobarconn(10) = sum(junk1(:))./ length(find(junk1));
    
%% Combine the values for all five networks;
fcvalues = [lpbafcvalues,lobarconn];

%% fill in the output variables;
fcvalues_all(jj,:) = fcvalues;

end

%% make the plot for the avg adjacency matrix
fcmatrix_avg = squeeze(mean(fcmatrix_all,1));
[a b] = sort(labelnums);
fcmatrix_trans = fcmatrix_avg(b,b);
% load the atlas labels
load('LPBA40labels.mat');
labels = labels(1:nroi);
labels_trans = labels(b);
f = figure; 
imagesc(fcmatrix_trans);
f.Position = [100 100 900 800]; 

ticklabels_new = cell(size(labels_trans));
for i = 1:length(ticklabels_new)
    if i<=14
        ticklabels_new{i} = ['\color{blue} ' strrep(labels_trans{i},'_',' ')];
    elseif i>14 & i<=26
        ticklabels_new{i} = ['\color{green} ' strrep(labels_trans{i},'_',' ')];
    elseif i>26 & i<=38
        ticklabels_new{i} = ['\color{orange} ' strrep(labels_trans{i},'_',' ')];
    else
        ticklabels_new{i} = ['\color{red} ' strrep(labels_trans{i},'_',' ')];
    end
end
set(gca, 'XTick',1:nroi, 'XTickLabel',ticklabels_new,'XTickLabelRotation',90);    
set(gca, 'YTick',1:nroi, 'YTickLabel',ticklabels_new);    

colorbar; 
%caxis([0 0.7]);
switch foi
    case {'gamma'}
        caxis([0.1 0.8]);
    otherwise
        caxis([0.1 0.8]);
end

figname = ['Average Age ' num2str(age) ' fcanalysis_' methodtype '_' foi '_' parcmethod '_' gridresolution '_' fcmethod '_' atlastype '.bmp'];
saveas(f,[ouputfolder figname]);
% colorbar('Ticks',[-5,-2,1,4,7],...
%          'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})

%% Output the data into a text file;
% 1) output the average matrix file
textfile     = ['Average Age ' num2str(age) ' fcanalysis_' methodtype '_' foi '_' parcmethod '_' gridresolution '_' fcmethod '_' atlastype '.edge'];
outputname   = [ouputfolder textfile];
fcmatrix_avg = squeeze(mean(fcmatrix_all,1));
dlmwrite(outputname,fcmatrix_avg);
% 2) output the stats file by lobar regions
textfile    = ['All Participants Age ' num2str(age) ' fcanalysis_' methodtype '_' foi '_' parcmethod '_' gridresolution '_' fcmethod '_' atlastype '.txt'];
outputname  = [ouputfolder textfile];
if exist(outputname) & replacearg == 0
    disp([outputname ' already exists. Delete it or replace it.'])
    return
end

fband_cells = cell(nsubjs,1);
fband_cells(:) = {foi};
subjects = participantnumbers.(['age' num2str(age)])';
subjages = ones(nsubjs,1)*age;
T=[];
T=table(subjects,subjages,fband_cells,fcvalues_all(:,1),fcvalues_all(:,2),fcvalues_all(:,3),fcvalues_all(:,4),fcvalues_all(:,5),fcvalues_all(:,6),fcvalues_all(:,7),fcvalues_all(:,8),fcvalues_all(:,9),fcvalues_all(:,10),fcvalues_all(:,11)...
,'VariableNames',{'participantnumber','age','freqband','lpba','FF','FT','FP','FO','TT','TP','TO','PP','PO','OO'});
size(T)
writetable(T,outputname,'Delimiter','\t');


