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

%% Set parameters and paths
participantnumber = 2;
cfg     = [];
cfg.age = 36;
cfg.eegfilepath = fullfile(folderBase,'Data',['Age' num2str(cfg.age) 'mos/']);
cfg.filepath    = fullfile(folderBase,'Outputs/'); % output file path;
cfg.modelspath  = fullfile(folderBase,'Sourcemodels/');
cfg.modelspath  = '/Users/wanzexie/Dropbox/Collaboration/SourceSpaceFC/Sourcemodels/';
% define the name of the eeg input file
eegfilename = ['Experiment 1 Subject ' num2str(participantnumber) ' Age ' num2str(cfg.age) ' Fieldtrip 1Hz Highpass_New_2s.mat']; 
cfg.eegfilename = eegfilename;
% The example data that you can download from my Dropbox folder are in (.mat) Fieldtrip preprocessing format.
% The SourceSpaceFC_PPS program should also work with data in other format, as long as it can be read into Fieldtrip.

%% Main
% phase-to-phase synchrony
fois = {'theta','alpha','beta','gamma'};
for i = 1:length(fois)
    cfg.foi = fois{i};
    % phase-to-phase synchrony
    SourceSpaceFCanalysis_PPS(cfg, participantnumber);
    % orthogonalized amplitude-to-amplitude correlation
    SourceSpaceFCanalysis_AAC(cfg, participantnumber);
end