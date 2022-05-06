function [icarejcomps EEG_ICA_edited] = DefineArtificialICAs_AdjustAndSASICA(cfg,EEG)
% This program aims to isolate artificial ICAs, e.g., VEOG, HEOG, and focal component, from the EEG data.
% However, "No automated method can accurately isolate artifacts without supervision (Chaumon et al., 2015),"
% which I totally agree; Therefore, if you don't trust the following automatic pipeline please use SASICA or 
% Adjust or other toolboxes to visually inspect all the ICAs, e.g., "PlotAndCheckICAsWithAdjust.m";
if ~isfield(cfg, 'focalcomp'),        cfg.focalcomp.enable = 1;                       end
if ~isfield(cfg, 'EOGcorr'),          cfg.EOGcorr.enable   = 1;                       end
if ~isfield(cfg, 'autocorr'),         cfg.autocorr.enable  = 1;                       end
if ~isfield(cfg, 'opts'),             cfg.opts.noplot      = 1;                       end
if ~isfield(cfg, 'plotarg'),          cfg.plotarg          = 0;                       end; plotarg = cfg.plotarg;          
%% Load the EEG data with ICA weights;
%EEG1 = pop_eegchanoperator( EEG, { 'ch125=(ch21+ch14+ch15-ch12-ch5-ch6)/3 label VEOG' 'ch126=(ch25+ch32+ch26-ch8-ch2-ch1)/3 label HEOG' }); 

%% run the ICAs detection with SASICA and Adjust with the program modified by WX;
% EOG correlation test;
cfg.EOGcorr.corthreshH    = 'auto 2.5'; % M +/- 2.5*SD (99%); M +/- 2*SD (95%);
cfg.EOGcorr.Veogchannames = 'VEOG';
cfg.EOGcorr.corthreshV    = 'auto 2.5';
cfg.EOGcorr.Heogchannames = 'HEOG';

% Focal component test;
cfg.focalcomp.focalICAout  = 'auto 3'; 

% Muscle movement (auto  correlation);
cfg.autocorr.dropautocorr = 'auto 3';

% ADJUST
cfg.ADJUST.enable = 1;
%[EEG cfg] = WX_eeg_SASICA(EEG,cfg);
[EEG1 output] = WX_eeg_SASICA(EEG,cfg);
%% Remove the bad components;
% Results from SASICA
%focalcomps  = find(EEG.reject.SASICA.icarejfocalcomp); 
%eogcomps    = find(EEG.reject.SASICA.icarejchancorr);
focalcomps  = find(EEG1.reject.SASICA.icarejfocalcomp); 
eogcomps    = find(EEG1.reject.SASICA.icarejchancorr);

if  cfg.autocorr.enable  == 0;
    mclcomps = [];
else
    mclcomps    = find(EEG1.reject.SASICA.icarejautocorr);
    gdsf_adjust = find(EEG1.reject.SASICA.icaADJUST.GDSF > EEG1.reject.SASICA.icaADJUST.soglia_GDSF);
    mclcomps    = intersect(mclcomps,gdsf_adjust);
end
% all bad components in SASICA
rejcomps_sasica = union(focalcomps, eogcomps);

% V Eye movements from ADJSUT:
% 1. Maximum Epoch Variance (MEV)
% 2. Spatial Average Difference (SAD)
mevcomps = find(EEG1.reject.SASICA.icaADJUST.nuovaV > EEG1.reject.SASICA.icaADJUST.soglia_V);
sadcomps = find(EEG1.reject.SASICA.icaADJUST.SAD > EEG1.reject.SASICA.icaADJUST.soglia_SAD);
vecomps_adjust = intersect(mevcomps,sadcomps);
% H Eye movements from ADJSUT:
% 1. MEV 2.Spatial Eye Difference (SED);
sedcomps = find(EEG1.reject.SASICA.icaADJUST.SED > EEG1.reject.SASICA.icaADJUST.soglia_SED);
hecomps_adjust = intersect(mevcomps,sedcomps);
% Blink
% 1. SAD 2. Spatial Variance Difference (SVD)
blkcomps = EEG1.reject.SASICA.icaADJUST.blink; % This blink detection with infant data in ADJSUT is really shitty; 
% over all eye related components by ADJUST;
eyecomps_adjust = unique([vecomps_adjust,hecomps_adjust,blkcomps]);
% Focal (Discontinuities) components from ADJUST;
gdsfcomps_adjust = EEG1.reject.SASICA.icaADJUST.disc;
% find out the overlapping problematic componencts by sasica and ADJUST;
icarejcomps_eye   = intersect(eogcomps,eyecomps_adjust);
icarejcomps_focal = intersect(focalcomps,gdsfcomps_adjust); % The focal comps detection is faily consistent between SASICA and ADJSUT;
icarejcomps       = union(icarejcomps_eye,icarejcomps_focal);
% add the mucle movement components to the icarejcomps array;
if cfg.autocorr.enable == 1 % it means auto correlation has been run and we get a 1 x 124 matrix for the results;
    icarejcomps = union(icarejcomps,mclcomps);
end

% Both SASICA and ADJUST save their results in the gcompreject field, which is the default place for bad comps in EEGLAB;
EEG1.reject.gcompreject=EEG1.reject.gcompreject.*0;  % clear ADJUST's results;           
EEG1.reject.gcompreject(icarejcomps) = 1;  % Substitute it with the final rejected comps for the dataset;
disp(['Remove ' num2str(length(icarejcomps)) ' components: ' num2str(icarejcomps)]);
EEG_ICA_edited = pop_subcomp( EEG1, icarejcomps, 0);
if plotarg;
   % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','on'); 
   % eeglab redraw;
    pop_eegplot(EEG1,1,0,0,'','dispchans',10,'spacing',100);
    set (gcf,'name','Plot of EEG data before ICA rejection');
    EEG_ICA_edited.reject = EEG1.reject;
    pop_eegplot(EEG_ICA_edited,1,0,0,'','dispchans',10,'spacing',100,'children',gcf);
    set (gcf,'name','Plot of EEG data after ICA rejection');
    disp(['These components are removed: ' num2str(icarejcomps)]);
    SASICA(EEG);
end
clear EEG1 


