function [] = SourceSpaceFCanalysis_AAC(cfg, participantnumber);
% This program will run source-space FC analysis using orthoganolized 
% amplitude to amplitude correlation (AAC).
% It is modified based on SourceSpaceFCanalysis_PPS, and please type
% "help SourceSpaceFCanalysis_PPS" for the description of the parameters.

%% Define parameters and methods;                   
if ~isfield(cfg, 'foi'),              cfg.foi              = 'alpha';                 end; foi        = cfg.foi;
if ~isfield(cfg, 'fcmethod'),         cfg.fcmethod         = 'AAC';                   end; fcmethod   = cfg.fcmethod;
if ~isfield(cfg, 'atlastype'),        cfg.atlastype        = 'LPBA';                  end; atlastype  = cfg.atlastype;
if ~isfield(cfg, 'replace'),          cfg.replacearg       = 0;                       end; replacearg = cfg.replacearg;
if ~isfield(cfg, 'parcmethod'),       cfg.parcmethod       = 'average';               end; parcmethod = cfg.parcmethod; 
if ~isfield(cfg, 'methodtype'),       cfg.methodtype       = 'eloreta';               end; methodtype = cfg.methodtype;
if ~isfield(cfg, 'gridresolution'),   cfg.gridresolution   = '6mm';                   end; gridresolution  = cfg.gridresolution;
if ~isfield(cfg, 'plotarg'),          cfg.plotarg          = 0;                       end; plotarg = cfg.plotarg;
cfg_original = cfg;
%% Global variables;
global EEG_FT  sourcedata moments m roitrialdata  
global filter atlas roivol
ft_defaults;

%% Define file names and Load the EEG data
age = cfg.age;
if age > 12 %
    segmentedtype='segmented'; % FEM model;
else
    segmentedtype='segmented_nma';% FEM model but with non-myelinated axons;
end
switch age   
    case 12
        mrinumber = 306; 
    case 36
        mrinumber = 328;
end
if mrinumber>=1000
    mristring=num2str(mrinumber);
elseif mrinumber>100
    mristring = ['0' num2str(mrinumber)];
else 
    mristring =  ['00' num2str(mrinumber)];
end

% define paths
eegfilepath = cfg.eegfilepath;
filepath    = cfg.filepath;
modelspath  = cfg.modelspath;

% check if the output matrix already exists;
outputname = ['Participant ' num2str(participantnumber) ' Age ' num2str(age) ' fcanalysis_' methodtype '_' foi '_' parcmethod '_' gridresolution '_' fcmethod '_' atlastype '.mat'];
if exist([filepath outputname]) & replacearg == 0;
    disp(['The following file already exists:' outputname])
    return
end

%load the EEG data; 
eegfilename = cfg.eegfilename;
load ([eegfilepath eegfilename],'EEG_FT'); 
%change the FT format from preprocessing to timelock format, which makes some of the following analyses easier;
cfg = [];
cfg.keeptrials = 'yes';
cfg.covariance = 'yes';
EEG_FT = ft_timelockanalysis(cfg,EEG_FT);

%% desampling to 250 Hz:
disp('Resampling the data to 250 Hz ...');
cfg=[];
cfg.resamplefs = 250;
EEG_FT = ft_resampledata(cfg, EEG_FT);

%% concatenate the data
catdata = permute(EEG_FT.trial,[2 3 1]); % change the dim to chan x samples x trials;
catdata = reshape(catdata,size(catdata,1),[],1); % concatenate all trials;

%% Filter the data before source analyis, suggested by Jessica Greeen and used in Tokariev et al. 2019.
% filter the data using matlab "bandpass" function in the signal processing toolbox;
disp(['Filter the data into different frequency bands using the Matlab "bandpass" function']);
fs = 250;
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
[catdata_filt d]= bandpass(catdata',bpfreq,fs,'ImpulseResponse','fir','Steepness',0.85);%'Steepness',0.85 is the default;
% disp the filter and output;
if plotarg
    fvtool(d);
    figure;
    FLIMS = [0 40];
    pspectrum([catdata_filt],fs,'FrequencyLimits',FLIMS);
    legend(['Steepness = 0.85;bp = [' num2str(bpfreq(1)) ' ' num2str(bpfreq(2)) ']'],'Location','south');
    xline(bpfreq(1));xline(bpfreq(2));
end

catdata_filt = catdata_filt';

%% Cortical source reconstruction;
mrifieldtripfoldername =  modelspath;
% If use eLORETA then there is no need to do source analysis from the scratch. Instead, use the filter created for each average MRI;
% However, if use beamforming or MNE then the source analysis needs to be conducted for each dataset because
% their source analysis calculates the covariance matrix based on each dataset;
% The name of the electrodetype needs to be defined regardless of the methodtype for source analysis;
electrodetype='hydrocelgsn124';
%load the grid; %The grid will be used later for atlas stuff regardless of source analysis methods;
filename=['S' mristring '_sourcemodel_' gridresolution '_ft_grid.mat'];
if ~exist([mrifieldtripfoldername filename]);
 disp(['does not exist quick return ' mrifieldtripfoldername filename]);
 return
end
disp(['load ' mrifieldtripfoldername filename]);
load([mrifieldtripfoldername filename],'grid');
 
filedirectory_lf = modelspath;
filename_filter=['S' mristring '_' electrodetype '_' segmentedtype '_' gridresolution '_ft_' methodtype '_filter.mat'];

if exist([filedirectory_lf filename_filter]);
    disp(['load ' filedirectory_lf filename_filter]);
    load([filedirectory_lf filename_filter],'filter');
    efilter = filter;
elseif ~exist([filedirectory_lf filename_filter]) & strcmp(methodtype, 'MNE');
    disp(['Create MNE filter for participant ' num2str(participantnumber) '...........']);
    EEG_FT_MNE = EEG_FT;
    [ntrial nchan nsamp] = size(EEG_FT.trial);
    segdata = reshape(catdata_filt,nchan,ntrial,nsamp); % segment the data into trials;
    segdata = permute(segdata,[2 1 3]); % change the dim to trial x chan x samples;
    EEG_FT_MNE.trial = segdata;
    
    cfg = [];
    cfg.covariance = 'yes';
    tlckFC = ft_timelockanalysis(cfg, EEG_FT);
    labels = tlckFC.label;
    for i = 1:length(labels)
        label = labels{i};
        label(1) = []; %delete 'E'
        labels{i}= label;
    end
    tlckFC.label = labels;
    
    % source analysis
    nelec  = 124;
    lfname = ['S0' num2str(mrinumber) '_hydrocelgsn' num2str(nelec) '_' segmentedtype '_' num2str(gridresolution) '_ft_lf.mat'];
    volname = ['S0' num2str(mrinumber) '_' segmentedtype '_ft_vol.mat'];
    if ~exist([modelspath lfname]) | ~exist([modelspath volname])
        disp([lfname ' or ' volname  ' does not exist']);
        return
    end
    load([modelspath lfname],'lf');
    % load the headmodel 
    load([modelspath volname],'vol');
    % Load the grid
    cfg               = [];
    cfg.method        = 'mne';
    % cfg.grid          = leadfield;
    % cfg.headmodel     = headmodel;
    cfg.grid          = lf;
    cfg.headmodel     = vol;
    cfg.keepfilter    = 'yes'
    cfg.mne.prewhiten = 'yes';
    cfg.mne.lambda    = 3;
    cfg.mne.scalesourcecov = 'yes';
    sourceFC          = ft_sourceanalysis(cfg,tlckFC);
    efilter = sourceFC.avg.filter;
    disp(['MNE filter created for participant ' num2str(participantnumber) '...........']);
    clear lf vol segdata EEG_FT_MNE;
else 
    disp(['no saved eLORETA filter, quick return ' filedirectory_lf filename_filter]);
    return
end

sourcedata = [];
ntrial  = size(EEG_FT.trial,1);
ndipole = find(~cellfun(@isempty,efilter)); %This array equals find(grid.inside);
%ndipole = find(grid.inside);
ntime   = size(EEG_FT.time,2);
% sourcedata.trial = NaN(ntrial,length(ndipole),ntime);
sourcedata.trial = NaN(length(ndipole),ntime*ntrial);
sourcedata1.trial = NaN(length(ndipole),ntime*ntrial);

tic;
disp('Source analysis starts...');
for j = 1:length(ndipole)
    % fprintf('%.0f out of %.0f \n',j,length(ndipole));
    % moments = efilter{ndipole(j)} * catdata_filt; % catdata_filt: nchan x all samples;
    moments = efilter{ndipole(j)} * catdata_filt;
    [u, s, v] = svd(moments, 'econ');  %Find the dimension explaining/representing the most variance; 
    m = u(:,1)'* moments; 
    sourcedata.trial(j,:) = m;    
    if (mod(j,1000) == 1)
        textprogressbar(j);
    end
end
analysistime = toc;
disp(['The Source analysis took ' num2str(analysistime) 's']);
savesourcedata = 0;
if savesourcedata % these sourcedata are big in size; 
    save([sourcepath sourcefilename],'sourcedata');
else
end

%% Hilbert transformation
tic
disp('Hilbert transformation ...');
sourcedata.trial = hilbert(sourcedata.trial')';
hilberttime = toc;
disp(['The hilbert transformation took ' num2str(hilberttime) 's']);

%% Parcellation into ROIs;
tic;
%load the atlas;
switch atlastype
    case 'LPBA'
        atlasname = ['S' mristring '_LPBA40atlas_ft.mat'];
        disp(['load ' modelspath atlasname]);
        load([modelspath atlasname],'atlas');
        nroi = 56;
    case 'YeoFiveNetworks' % I am working on creating the age-ralted version of the Yeo Network atlas;
        YeoFiveNetworks = ['S0' num2str(mrinumber) '_YeoFiveNetworks_ft_atlas.mat'];
        disp(['load ' atlasespath YeoFiveNetworks]);
        load([atlasespath YeoFiveNetworks],'atlas');
        nroi = 30;
end

%load the atlas/brain ac info
brain_ac = ['S' mristring '_brain_atlas_ac.mat'];
load([modelspath brain_ac],'brain_ac');
%Find the dipoles;
positioninside=grid.pos(find(grid.inside),:);
positioninside=positioninside + repmat(brain_ac,size(positioninside,1),1);
indicesinside=sub2ind(atlas.dim,positioninside(:,1),positioninside(:,2),positioninside(:,3));

%Create representative data for each brain ROI;
for j=1:length(atlas.labelnumbers)
    % not use the subcortical regions;
    if strcmp(atlastype, 'LPBA') & j>48;
     continue
    end
    searchloc=find(atlas.inside == atlas.labelnumbers(j));    %indices for this label/ROI
    searchpos=atlas.pos(searchloc,:);                 %locations for this label/ROI
    searchpos=double(int64(searchpos))+repmat(brain_ac,size(searchpos,1),1);    
    searchindices=sub2ind(atlas.dim,searchpos(:,1),searchpos(:,2),searchpos(:,3));
    if strcmp(parcmethod,'centroid')
     % find out the centroid position;
     centroidpos = round(mean(searchpos));
     % create a centroid cube;
     X = [centroidpos(1)-10:centroidpos(1)+10];
     Y = [centroidpos(2)-10:centroidpos(2)+10];
     Z = [centroidpos(3)-10:centroidpos(3)+10];
     cubepos = combvec(X,Y,Z)';
     cubeindices=sub2ind(atlas.dim,cubepos(:,1),cubepos(:,2),cubepos(:,3));
     %get an array that marks every point in sourcedata with atlas (cubes)
     [jj kk]=ismember(indicesinside,cubeindices);
     [junk1 junk2]=ismember(indicesinside,searchindices);
     junk3 = sum([jj junk1],2);
     % junk3 == 2;
     jj(junk3~=2) = 0;
     disp([atlas.labels{j} ' has ' num2str(sum(jj)) ' centroid source voxels.']);
     nofvoxels(j) = sum(jj);
    else
      %get an array that marks every point in sourcedata with atlas (ROIs)
     [jj kk]=ismember(indicesinside,searchindices);
     disp([atlas.labels{j} ' has ' num2str(sum(jj)) ' source voxels.']);
    end
    atlasmatrix(j,:)=jj';  
end
roivol=atlasmatrix*repmat(1,[size(atlasmatrix,2) 1]);    %sum over columns
atlasmatrixdiv=atlasmatrix./repmat(roivol,[1 size(atlasmatrix,2)]); % This is for the average method;

roitrialdata=[];
if strcmp(atlastype, 'LPBA')
    roitrialdata.trial = NaN(48,ntime*ntrial);
else
    roitrialdata.trial = NaN(length(atlas.labelnumbers),ntime*ntrial);
end
switch parcmethod
    case 'PCA'
        for j = 1:length(atlas.labelnumbers)
            % not use the subcortical regions;
            if strcmp(atlastype, 'LPBA') & j>48;
                continue
            end
            roidata = [];
            roidata = sourcedata.trial(find(atlasmatrix(j,:)),:);
            roidata = squeeze(roidata)';
            [coeff score latent useless explained] = pca(roidata);
            roidata_maxpca =  score(:,1)'; % pick up the pca component explaining the most of the variance;
            %{ 
            roidata_pcas = roidata*coeff;  %original component, and the scores are these components minus the means;
            roidata_maxpca2 = roidata_pcas(:,1)'; 
            %} 
            roitrialdata.trial(j,:)=roidata_maxpca;
        end
    case {'average','centroid'}
        roitrialdata.trial(:,:) = atlasmatrixdiv*sourcedata.trial(:,:);
end

sample = size(roitrialdata.trial,2);
roitrialdata.time = [0.001:(1/250):(sample./250)];
roitrialdata.dimord = 'chan_time';
roitrialdata.label = atlas.labels;

analysistime = toc;
disp(['The averaging/PCA/centroid process took ' num2str(analysistime) 's']);

%% FC analysis with orthogonalized amplitude & amplitude correlation;
[Rp, Ro] = OrthogonalPowCorr(roitrialdata.trial);
if plotarg
    figure; imagesc(Rp); colorbar; caxis([0 0.4]);
    figure; imagesc(Ro); colorbar; caxis([0 0.2]);
end

%% save source fc matrices
disp('Save source data & FC matrices ...');
fcmatrix.Rp = Rp;
fcmatrix.Ro = Ro;
save([filepath outputname],'fcmatrix');

% you may also save the original cfg structure with all the file names, paths, and parameters;
% comment out or delete the following two lines if you don't want to save the cfg;
outputname_cfg = strrep(outputname,'.mat','_cfg.mat');
save([filepath outputname_cfg],'cfg_original');
