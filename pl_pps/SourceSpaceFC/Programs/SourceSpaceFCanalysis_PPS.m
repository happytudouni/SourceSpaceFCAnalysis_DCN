function [] = SourceSpaceFCanalysis_PPS(cfg, participantnumber)
% run this program under the "Programs/" directory OR change the paths to local directories accordingly;

%% Define parameters and methods;                   
if ~isfield(cfg, 'foi'),              cfg.foi              = 'alpha';                 end; foi        = cfg.foi;
if ~isfield(cfg, 'fcmethod'),         cfg.fcmethod         = 'wpli';                  end; fcmethod   = cfg.fcmethod;
if ~isfield(cfg, 'atlastype'),        cfg.atlastype        = 'LPBA';                  end; atlastype  = cfg.atlastype;
if ~isfield(cfg, 'replacearg'),       cfg.replacearg       = 0;                       end; replacearg = cfg.replacearg;
if ~isfield(cfg, 'parcmethod'),       cfg.parcmethod       = 'centroid';              end; parcmethod = cfg.parcmethod; 
if ~isfield(cfg, 'methodtype'),       cfg.methodtype       = 'eloreta';               end; methodtype = cfg.methodtype;
if ~isfield(cfg, 'gridresolution'),   cfg.gridresolution   = '6mm';                   end; gridresolution  = cfg.gridresolution;
if ~isfield(cfg, 'hztype'),           cfg.hztype           = 1;                       end; hztype  = cfg.hztype;
if ~isfield(cfg, 'age'),              cfg.age              = 36;                      end; age     = cfg.age;
if ~isfield(cfg, 'plotarg'),          cfg.plotarg          = 0;                       end; plotarg = cfg.plotarg;

%% Global variables;
global EEG_FT  sourcedata moments m roitrialdata  
global filter atlas roivol
ft_defaults;

%% Define file names and Load the EEG data
% The source models are slightly different for children under and above 12 months of age;
% Specially, for children <= 12 mos their FEM model includes non-myelinated axons, whereas for children > 12 mos the FEM model does not include non-myelinated axons; 
if age > 12 
    segmentedtype='segmented'; % FEM model;
else
    segmentedtype='segmented_nma';% FEM model but with non-myelinated axons;
end
% The following mrinumbers represent mri templates (downloaded from the JERLAB at U of South Carolina) for children at different ages ;
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
filepath = '../Outputs/';
modelspath  = '../Sourcemodels/';
eegfilepath =  ['../Data/Age' num2str(age) 'mos/'];
sourcepath  = '../SourceData/';

% check if the output matrix already exists;
outputname = ['Participant ' num2str(participantnumber) ' Age ' num2str(age) ' fcanalysis_' methodtype '_' foi '_' parcmethod '_' gridresolution '_' fcmethod '_' atlastype '.mat'];
if exist([filepath outputname]) & replacearg == 0;
    disp(['The following file already exists:' outputname])
    return
end

%load the EEG data; 
eegfilename = ['Experiment 1 Subject ' num2str(participantnumber) ' Age ' num2str(age) ' Fieldtrip ' num2str(hztype) 'Hz Highpass_New_2s.mat']; 
load ([eegfilepath eegfilename],'EEG_FT'); 
% change the data format from preprocessing to timelock format, which makes some of the following analyses easier;
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

%% filter the data using matlab "bandpass" function in the signal processing toolbox;
disp(['Filter the data into different frequency bands using the Matlab "bandpass" function']);
fs = 250; % sampling rate
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
% However, if use beamforming then the source analysis needs to be conducted for each dataset because
% beamforming source analysis calculates the covariance matrix for each dataset.;
% We provided the eLORETA filter along with this program and we recommend using the eLORETA filter unless you want to try to create a Beamformer filter in Fieldtrip;
% Instructions on source analysis are out of the scope of the current tutorial. Please see tutorials on the Fieldtrip website, e.g.,: https://www.fieldtriptoolbox.org/tutorial/beamformer/;

% There were 124 (not 128) electrodes for infants and young children's EGI nets (HGSN);
electrodetype='hydrocelgsn124';

%load the grid; %The grid will be used later for atlas stuff regardless of source analysis methods;
filename=['S' mristring '_sourcemodel_' gridresolution '_ft_grid.mat'];
if ~exist([mrifieldtripfoldername filename]);
 disp(['does not exist quick return ' mrifieldtripfoldername filename]);
 return
end
 disp(['load ' mrifieldtripfoldername filename]);
 load([mrifieldtripfoldername filename],'grid');
 
% load the eLORETA filter 
filedirectory_lf = modelspath;
filename_filter=['S' mristring '_' electrodetype '_' segmentedtype '_' gridresolution '_ft_' methodtype '_filter.mat']
if exist([filedirectory_lf filename_filter]);
    disp(['load ' filedirectory_lf filename_filter]);
    load([filedirectory_lf filename_filter],'filter');
    efilter = filter;
else
    disp(['no saved filter quick return ' filedirectory_lf filename_filter]);
    return
end

% create an empty source data structure
sourcedata = [];
ntrial  = size(EEG_FT.trial,1);
% find the index of the dipoles in the efilter structure because there are some empty cells in the efilter structure;
ndipole = find(~cellfun(@isempty,efilter)); %This array equals find(grid.inside);
% figure out the length of each epoch;
ntime   = size(EEG_FT.time,2);
% create the empty sourcedata structure with nans. This name of this structure ("sourcedata") is arbitrary;
sourcedata.trial = NaN(length(ndipole),ntime*ntrial);

% The start of the source reconstructure
tic;
disp('Source analysis starts...');
textprogressbar('initialize');
textprogressbar(['num of dipoles finished  ',length(ndipole)]);
for j = 1:length(ndipole)
    % do the source reconstruction dipole by dipole
    % moments represents the source activation in one dipole, and it has 3 orientations/dimensions;
    moments = efilter{ndipole(j)} * catdata_filt;
    % Use Singular Value Decomposition to find out the virtual dimension that explains the most variance of the dipole activation; 
    [u, s, v] = svd(moments, 'econ');  
    m = u(:,1)'* moments; 
    % m represents the "strongest" dimension or orientation of the dipole activation, so the 3-D dipole activation is reduced to 1-D;
    sourcedata.trial(j,:) = m;    
    if (mod(j,1000) == 1)
        textprogressbar(j);
    end
end
analysistime = toc;
textprogressbar('finished');
disp(['The Source analysis took ' num2str(analysistime) 's']);

%% Parcellation of the dipoles across the brain into ROIs;
tic;
%load the brain atlas;
switch atlastype
    case 'LPBA'
        atlasname = ['S' mristring '_LPBA40atlas_ft.mat'];
        disp(['load ' modelspath atlasname]);
        load([modelspath atlasname],'atlas');
        nroi = 56;
end

%load the atlas/brain anterior commissure info
brain_ac = ['S' mristring '_brain_atlas_ac.mat'];
load([modelspath brain_ac],'brain_ac');

% Find the postions of the dipoles with activity;
positioninside=grid.pos(find(grid.inside),:);
% add the 3D coordinates of the AC to the positions of the dipoles to align the dipoles with the brain atlas;
positioninside=positioninside+double(int64(repmat(brain_ac,size(positioninside,1),1)));
indicesinside=sub2ind(atlas.dim,positioninside(:,1),positioninside(:,2),positioninside(:,3));

%Create representative data for each brain ROI;
% from Line 199 to line 256 are the codes for parcellation of dipoles, i.e., calculate the source activation for each ROI in the LPBA atlas;
for j=1:length(atlas.labelnumbers)
 if strcmp(atlastype, 'LPBA') & j>48; % in the LPBA atlas, subcortical regions have label numbers > 48, so skip these subcortical regions;
     continue
 end
 searchloc=find(atlas.inside == atlas.labelnumbers(j));    %indices for this label/ROI
 searchpos=atlas.pos(searchloc,:);                 %locations for this label/ROI
 searchpos=double(int64(searchpos))+repmat(brain_ac,size(searchpos,1),1);    
 searchindices=sub2ind(atlas.dim,searchpos(:,1),searchpos(:,2),searchpos(:,3));
 % Use the centroid method for parcellation
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
     disp([atlas.labels{j} ' has ' num2str(sum(jj)) ' centroid voxels.']);
 else
      %get an array that marks every point in sourcedata with atlas (ROIs)
     [jj kk]=ismember(indicesinside,searchindices);
 end
 atlasmatrix(j,:)=jj';  
end
roivol=atlasmatrix*repmat(1,[size(atlasmatrix,2) 1]);    %sum over columns
atlasmatrixdiv=atlasmatrix./repmat(roivol,[1 size(atlasmatrix,2)]); 

roitrialdata=[];
if strcmp(atlastype, 'LPBA')
    roitrialdata.trial = NaN(48,ntime*ntrial);
    roitrialdata.label = atlas.labels(1:48);
else
    roitrialdata.trial = NaN(length(atlas.labelnumbers),ntime*ntrial);
end
switch parcmethod
    case 'PCA'
        for j = 1:length(atlas.labelnumbers)
            roidata = [];
            roidata = sourcedata.trial(find(atlasmatrix(j,:)),:);
            roidata = squeeze(roidata)';
            [coeff score latent useless explained] = pca(roidata);
            roidata_maxpca =  score(:,1)'; % pick up the pca component explaining the most of the variance;
            roitrialdata.trial(j,:)=roidata_maxpca;
        end
    case {'average','centroid'}
        roitrialdata.trial(:,:) = atlasmatrixdiv*sourcedata.trial(:,:);
end

sample = size(roitrialdata.trial,2);
roitrialdata.time = [0.001:(1/250):(sample./250)];
roitrialdata.dimord = 'chan_time';

analysistime = toc;
disp(['The averaging/PCA/centroid process took ' num2str(analysistime) 's']);

%% segment the concatenated source data back into epochs
% We are doing this because wPLI and other PPS methods prefer data in trials/epochs;
% Please see the following tutorial for FC analysis in Fieldtrip: https://www.fieldtriptoolbox.org/tutorial/connectivity/
data     = roitrialdata.trial;
sampling = 250;
nch      = size(data,1);
ntrials  = size(data,2)./(250*2); % 2s epoch
data     = reshape(data,[nch sampling*2 ntrials]); % chan_time_trial(rpt)
data     = permute(data,[3 1 2]); % change the dim to 'rpt_chan_time';
roitrialdata.trial  = data;
roitrialdata.dimord = 'rpt_chan_time';
roitrialdata.time   = EEG_FT.time;
roitrialdata.sampleinfo = EEG_FT.sampleinfo;
nroi     = size(roitrialdata.trial,2);

%% repetition / permutation
permutation = 1; % change this to 0 if repetition is not needed;
if permutation & strcmp(fcmethod,'wpli') & size(roitrialdata.trial,1)>=30
    ntrial      = size(roitrialdata.trial,1); % total number of trials
    nrep        = 50; % number of permutation
    ntr_round   = 30; % number of trials selected for each round

    lowernumber   = ntr_round; % number of trials per round;
    biggernumber  = ntrial;
    % create a 50 (nrep) by 30 (ntr_round) matrix "a"
    % in this matrix, each row contains the idx of the trials (30) to be used for that round of FC calculation;
    a = zeros(nrep,lowernumber);
    for i = 1:nrep
        a(i,:) = sort(randperm(biggernumber,lowernumber)); 
    end
    fcmatrix    = []; 
    repfcmatrix = zeros(nrep,nroi,nroi);
    for i = 1:size(a,1);
        EEG_pm = roitrialdata;
        EEG_pm.trial = EEG_pm.trial(a(i,:),:,:);
    %% Frequency analysis
    % A few options for cfg.method in Fieldtrip: {'coh', 'csd', 'wpli', 'wpli_debiased', 'plv','imag'}
        cfg            = [];
        cfg.output     = 'fourier';
        cfg.method     = 'mtmfft';
        cfg.foilim     = bpfreq;
        cfg.keeptrials = 'yes';
        cfg.taper      = 'hanning';
        EEG_freq       = ft_freqanalysis(cfg, EEG_pm);

    %% FC analysis
        cfg = [];
        cfg.method = fcmethod;
        if strcmp(fcmethod,'imag');
            cfg.method = 'coh';
            cfg.complex = 'absimag';
        end
        source_conn = ft_connectivityanalysis(cfg, EEG_freq); % source_conn is in Fieldtrip output format

        wplispctrm  = source_conn.wplispctrm;
        wplispctrm  = abs(source_conn.wplispctrm); %use the abs value for wpli;
        wplispctrm  = mean(wplispctrm,3); % average across the frequency bins; skip this step if you want to plot FC spetrum density by bins;
        wplispctrm(logical(eye(size(wplispctrm))))   = 0; % change the diagonal values to 0s
        zwplispctrm = 0.5*log((1+wplispctrm)./(1-wplispctrm)); % fisher's z-transformation
        repfcmatrix(i,:,:) = zwplispctrm;
    end
    fcmatrix.wpli  = squeeze(mean(repfcmatrix,1));
    if plotarg; figure; imagesc(fcmatrix.wpli); colorbar; caxis([0 0.4]);end % plot the adjacency matrix 
else % no permutation for the other methods;
    fcmatrix  = [];
    % frequency analysis
    cfg            = [];
    cfg.output     = 'fourier';
    cfg.method     = 'mtmfft';
    cfg.foilim     = bpfreq;
    cfg.keeptrials = 'yes';
    cfg.taper      = 'hanning';
    EEG_freq       = ft_freqanalysis(cfg, roitrialdata);
    % connectivity analysis
    cfg = [];
    cfg.method = fcmethod;
    if strcmp(fcmethod,'imag');
        cfg.method = 'coh';
        cfg.complex = 'absimag';
    end
    source_conn = ft_connectivityanalysis(cfg, EEG_freq);
    switch fcmethod
        case {'coh','imag'} % please see ft_connectivityanalysis for more options
            cohspctrm = source_conn.cohspctrm;
            cohspctrm = mean(cohspctrm,3);
            cohspctrm(logical(eye(size(cohspctrm))))   = 0;
            zcohspctrm = 0.5*log((1+cohspctrm)./(1-cohspctrm));
            fcmatrix.coh = zcohspctrm; %calculate the mean fc for different frequency bins;
        case {'wpli'} 
            wplispctrm  = source_conn.wplispctrm;
            wplispctrm  = abs(source_conn.wplispctrm); %use the abs value for wpli;
            wplispctrm  = mean(wplispctrm,3); % average across the frequency bins; skip this step if you want to plot FC spetrum density by bins;
            wplispctrm(logical(eye(size(wplispctrm))))   = 0; % change the diagonal values to 0s
            zwplispctrm = 0.5*log((1+wplispctrm)./(1-wplispctrm)); % fisher's z-transformation
            fcmatrix.wpli = zwplispctrm; %calculate the mean fc for different frequency bins;
        otherwise
    end   
end  % if permutation & strcmp(fcmethod,'wpli')

disp('Save FC matrices ...');
save([filepath outputname],'fcmatrix');


