clear all
clc

%% Set environment
restoredefaultpath;

dir_matlab = fullfile('..','..','..','..','anwenkang@gmail.com - Google Drive','My Drive','Research_MATLAB/');
addpath(fullfile(dir_matlab,'eeglab2019_0/'));
addpath(fullfile(dir_matlab,'fieldtrip-master'));

eeglab;
close all;
ft_defaults;

%% Set parameters
cfg = [];
participantnumber = 2;

%% Set path
cfg.path_to_data = '';

%% Main
SourceSpaceFCanalysis_PPS(cfg, participantnumber);