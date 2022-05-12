
clear all
clc

%% Set environment
% set up your matlab path and do addpath for the eeglab and fieldtrip folders;
dir_matlab = fullfile('/','Users','wanzexie','Documents','MATLAB');
eeglabpath    = 'eeglab2021.1'; % define your eeglab path;
fieldtrippath = 'fieldtrip-20210409'; % define your fiedltrip path;
addpath(fullfile(dir_matlab,'toolbox',eeglabpath));
addpath(fullfile(dir_matlab,'toolbox',fieldtrippath));

eeglab;
close all;
ft_defaults;

% define the path of your SourceSpaceFC folder
folderBase  = fullfile('/','Users','wanzexie','Documents','GitHub','SourceSpaceFCAnalysis_DCN','SourceSpaceFC','pl_pps');
cd(folderBase);
programpath = fullfile(folderBase,'Programs/');
addpath(programpath);


folderBase  = fullfile('/','Users','wanzexie','Documents','GitHub','SourceSpaceFCAnalysis_DCN','SourceSpaceFC','pl_pps');
cd(folderBase);
programpath = fullfile(folderBase,'Programs/');
addpath(programpath);

%% Run analysis for all subjects;
ages = [12 36];
participantnumbers  = [];
for i = 1:length(ages)
    age = ages(i);
    datapath = fullfile(folderBase, 'Data', ['Age' num2str(age) 'mos/']);
    dirinfo = dir(datapath);    
    ii = 0;
    for K = 1 : length(dirinfo)
        if ~contains(dirinfo(K).name,'.mat')
            % matlab versions before R2016b does not have this "contains" function. Use 'strfind' instead;
        else
            ii = ii + 1;
            thisdir = dirinfo(K).name;
            thisstring = extractBetween(thisdir,"Subject "," Age");
            subjects(ii) = str2num(thisstring{1});
        end
    end
    participantnumbers.(['age' num2str(age)]) = unique(subjects);
end 

for i = 1:length(ages)
    age = ages(i);
    subjects = participantnumbers.(['age' num2str(age)]);
    for ii = 1:length(subjects)
    participantnumber = subjects(ii);
    % functional connectivity analysis  
    fois = {'theta','alpha','beta','gamma'};
    for iii = 1:length(fois)        
        cfg     = [];
        cfg.age = age;
        cfg.eegfilepath = fullfile(folderBase,'Data',['Age' num2str(cfg.age) 'mos/']);
        cfg.filepath    = fullfile(folderBase,'Outputs/'); % output file path;
        cfg.modelspath  = fullfile(folderBase,'Sourcemodels/');
        cfg.modelspath  = '/Users/wanzexie/Dropbox/Collaboration/SourceSpaceFC/Sourcemodels/';
        % define the name of the eeg input file
        eegfilename = ['Experiment 1 Subject ' num2str(participantnumber) ' Age ' num2str(cfg.age) ' Fieldtrip 1Hz Highpass_New_2s.mat']; 
        cfg.eegfilename = eegfilename;
        
        cfg.foi        = fois{iii};
        cfg.fcmethod   = 'wpli';    %'wpli','imag','coh' or other methods in ft_fcanalysis
        cfg.parcmethod = 'average'; %'centroid' or 'average'
        cfg.atlastype  = 'LPBA'; 
        cfg.methodtype = 'eloreta'; %'eloreta'  or 'MNE'
        cfg.replacearg = 0;
        cfg.plotarg    = 0;
        SourceSpaceFCanalysis_PPS(cfg, participantnumber);
    end
    end % for ii = 1:length(subjects)
end % for i = 1:length(ages)


disp(['The analysis has been conducted for Ages ' num2str(ages) '.']);