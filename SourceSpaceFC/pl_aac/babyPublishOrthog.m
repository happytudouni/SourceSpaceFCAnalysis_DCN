% This script describes the analysis procedures to produce frequency band-limited connectivity matrices
% from orthogonalized power envelopes derived from EEG. This script is annotated with comments throughout
% to aid in clarity.
% Conventions:
% functions prefixed by RTT_ are custom functions
% variables prefixed by thresh are values that are manually set specifically for these analyses
% fields in the EEG structure beginning with a capital letter are custom fields (ex. EEG.ChansBad is a custom field and
% EEG.data is a standard field)

%% 0. Preamble.

clearvars
close all
clc
% set up the eeglab default paths
eeglab;close;
%% 1. Generate unconstrained kernels in Brainstorm for each age group

% Brainstorm is an open-source toolbox for EEG and MEG available at https://neuroimage.usc.edu/brainstorm/.

% Brainstorm is used to generate the cortical source space model and compute the inverse solution, generating the
% imaging kernel. The following steps are applied for each age group:
% Refer to corresponding screenshots for a walk-through of the procedure.

% 1.1. Create protocol.
% 1.2. Define template anatomy.
% 1.3.a-b. Create cortex tesselation by downsampling the high-resolution template cortex.
% 1.4. Create a placeholder subject Plotter_5K to generate the imaging kernel from.
% 1.5.a-d. Set common files for the cap and noise modeling. Compute head model using OpenM/EEG plugin.
% 1.6.a-c. Import 1 second of EEG from any subject to serve as a template structure.
% 1.7.a-c. Compute the unconstrained imaging kernel for use in subsequent calculations.
% Save output as kernels.mat

%% 2. Calculate directional connectivity matrices from orthogonalized power envelopes

% user-defined base directory
folderBase = cd; % current directory

D = struct2table(dir(fullfile(folderBase,'PREPROC','FINAL')));
% D = D(~D.isdir,:);
D = D(contains(D.name,'set'),:);
D.PathAndName = fullfile(D.folder,D.name);

D.AgeSubject = extractBefore(D.name,'.set');
D.PathAndNameNewOPEC = strrep(strrep(D.PathAndName,'FINAL','OPEC'),'.set','.mat');

% load the unconstrained kernels for each age group
load(fullfile(folderBase,'ADMIN','masterKernelsInfantsToddlers5K.mat'),'KU*')
% load the frequency band delineations for each age group
load(fullfile(folderBase,'ADMIN','masterFrequencyBandsInfantsToddlers5K.mat'),'band*')
% declare the window duration coefficient
threshWindowCoeff = 10;

bands = {'THETA','ALPHA','BETA','GAMMA'};
for iBand=1:length(bands)
   band = bands{iBand};
   for iFile=1:height(D)
      varName = sprintf('RA_%s_BAND_%s',D.AgeSubject{iFile},band);
      if ~exist(D.PathAndNameNewOPEC{iFile},'file') || isempty(who('-file',D.PathAndNameNewOPEC{iFile},varName))
         % theta and alpha bands differ between the age groups
         if startsWith(D.AgeSubject{iFile},'A12')
            KU = KU12;
            bandTheta = bandThetaAge12;
            bandAlpha = bandAlphaAge12;
         elseif startsWith(D.AgeSubject{iFile},'A36')
            KU = KU36;
            bandTheta = bandThetaAge36;
            bandAlpha = bandAlphaAge36;
         else
            error('unknown age')
         end
         switch band
            case 'THETA'
               bandLimits = bandTheta;
            case 'ALPHA'
               bandLimits = bandAlpha;
            case 'BETA'
               bandLimits = bandBeta;
            case 'GAMMA'
               bandLimits = bandGamma;
            otherwise
               error('unknown band')
         end
         EEG = pop_loadset(D.PathAndName{iFile});
         cmw = RTT_CMW(bandLimits(1) + range(bandLimits) / 2,EEG.srate,range(bandLimits),'Hz',3);
         BPCD = RTT_Conv(EEG.data,cmw);
         % calculate the PCA reduced kernel
         KR = RTT_KernelPCA(KU,real(BPCD));
         clearvars KU
         [SSw,nWindows] = RTT_WindowData(KR * BPCD,cmw.length * threshWindowCoeff);
         RAw = single(NaN(size(SSw,1)));
         RAw = repmat(RAw,[1,1,nWindows]);
         for iWindow=1:nWindows
            RA = RTT_Orthog(squeeze(SSw(:,:,iWindow)));
            RA(~isfinite(RA)) = NaN;
            RAw(:,:,iWindow) = RA;
         end
         eval(sprintf('%s = nanmedian(RAw,3);',varName))
         eval(sprintf('%s = cmw;',regexprep(varName,'^RA_','cmw_')))
         if ~exist(D.PathAndNameNewOPEC{iFile},'file')
            eval(sprintf('save(''%s'',''*%s'')',D.PathAndNameNewOPEC{iFile},extractAfter(varName,'_')))
         else
            eval(sprintf('save(''%s'',''*%s'',''-append'')',D.PathAndNameNewOPEC{iFile},extractAfter(varName,'_')))
         end
         clearvars -except band* D KU* thresh* i* folderBase
      end
   end
end

disp('section 2 complete')

%% 3. Aggregate results

folderOPEC = fullfile(folderBase,'PREPROC','OPEC');

D = struct2table(dir(folderOPEC));
D = D(~D.isdir,:);
D.PathAndName = fullfile(D.folder,D.name);
D.name = extractBefore(D.name,'.');
D12 = D(startsWith(D.name,'A12'),:);
D36 = D(startsWith(D.name,'A36'),:);
% remove the file names with '_PRUNED' because they just .mat that contains the preprocessing ICA information;
D12(contains(D12.name,'PRUNED'),:) = [];
D36(contains(D36.name,'PRUNED'),:) = [];
clearvars D

bands = {'THETA','ALPHA','BETA','GAMMA'};
RAagg12 = single(NaN(5000,5000));
RAagg12 = repmat(RAagg12,[1,1,height(D12),length(bands)]);
RAagg36 = single(NaN(5000,5000));
RAagg36 = repmat(RAagg36,[1,1,height(D36),length(bands)]);

for iFile=1:height(D12)
   for iBand=1:length(bands)
      disp([12,iFile,iBand])
      band = bands{iBand};
      varName = sprintf('RA_%s_BAND_%s',D12.name{iFile},band);
      eval(sprintf('load(D12.PathAndName{iFile},''%s'');',varName))
      eval(sprintf('RAagg12(:,:,iFile,iBand) = %s;',varName))
      eval(sprintf('clearvars %s',varName))
   end
end

for iFile=1:height(D36)
   for iBand=1:length(bands)
      disp([36,iFile,iBand])
      band = bands{iBand};
      varName = sprintf('RA_%s_BAND_%s',D36.name{iFile},band);
      eval(sprintf('load(D36.PathAndName{iFile},''%s'');',varName))
      eval(sprintf('RAagg36(:,:,iFile,iBand) = %s;',varName))
      eval(sprintf('clearvars %s',varName))
   end
end

load(fullfile(folderBase,'ADMIN','masterCortexInfantsToddlers5K.mat'),'cortex*')
% thresh_oc is the underestimation correction factor as described in Hipp's 2012 Nature Neuroscience paper
thresh_oc = 0.57796;

[ROagg12 , RNagg12] = deal(RAagg12);
for iSubject=1:size(RAagg12,3)
   for iBand=1:size(RAagg12,4)
      disp([12,iSubject,iBand])
      
      RA = squeeze(RAagg12(:,:,iSubject,iBand));
      RO = RA ./ thresh_oc;
      RO = (RO + RO') ./ 2;
      ROagg12(:,:,iSubject,iBand) = RO;
      
      RN = (RO - median(RO,1,'omitnan')) ./ (1.4826 .* mad(RO,1,1));
      RNagg12(:,:,iSubject,iBand) = RN;
      
   end
end

[ROagg36 , RNagg36] = deal(RAagg36);
for iSubject=1:size(RAagg36,3)
   for iBand=1:size(RAagg36,4)
      disp([36,iSubject,iBand])
      
      RA = squeeze(RAagg36(:,:,iSubject,iBand));
      RO = RA ./ thresh_oc;
      RO = (RO + RO') ./ 2;
      ROagg36(:,:,iSubject,iBand) = RO;
      
      RN = (RO - median(RO,1,'omitnan')) ./ (1.4826 .* mad(RO,1,1));
      RNagg36(:,:,iSubject,iBand) = RN;
      
   end
end

disp('section 3 complete')

%% 4. Compute stats

Tagg12 = NaN(size(squeeze(ROagg12(:,:,1,:))));
for iBand=1:size(ROagg12,4)
   for iVertex=1:size(ROagg12,1)
      disp([12,iBand,iVertex])
      RNV = squeeze(RNagg12(:,iVertex,:,iBand));
      [t,p] = RTT_Tstat1(RNV,0,'b');
      Tagg12(:,iVertex,iBand) = t;
   end
end

Tagg36 = NaN(size(squeeze(ROagg36(:,:,1,:))));
for iBand=1:size(ROagg36,4)
   for iVertex=1:size(ROagg36,1)
      disp([36,iBand,iVertex])
      RNV = squeeze(RNagg36(:,iVertex,:,iBand));
      [t,p] = RTT_Tstat1(RNV,0,'b');
      Tagg36(:,iVertex,iBand) = t;
   end
end

finalResult12 = squeeze(mean(Tagg12,2));
finalResult36 = squeeze(mean(Tagg36,2));

save(fullfile(folderBase,'PREPROC','FINAL','finalResults.mat'),'finalResult*')

disp('section 4 complete')

%% custom functions

function cmw = RTT_CMW(cf,srate,bw,bwUnits,nSigmaT)
% constructs a complex Morlet wavelet for convolution
% cf = center frequency in Hz
% srate = sampling rate in samples/sec
% bw = total bandwidth in either Hz or octaves
% bwUnits = 'Hz' or 'octaves'
% nSigmatT = temporal std dev, reported as +/- (i.e., nSigmaT = +/- 3 == 6 total sigma)
   if startsWith(bwUnits,'h','IgnoreCase',true)
      % bandwidth is in units of Hz and describes whole range ([cf - bw/2,cf + bw/2])
      f1 = cf - bw/2;
      f2 = cf + bw/2;
      bwOctaves = log2(f2/f1);
   elseif startsWith(bwUnits,'o','IgnoreCase',true)
      % bandwidth is in units of octaves
      bwOctaves = bw;
      f1 = (2 * cf) / (2^bwOctaves + 1);
      f2 = 2^bwOctaves * f1;
   else
      error('bandwidth units must be either Hz or Octaves')
   end
   bwHz = f2 - f1;
   sigmaF = bwHz / 2;
   nCycles = cf / sigmaF;
   sigmaT = nCycles / (2 * pi * cf);
   t = 1/srate : 1/srate : nSigmaT*sigmaT;
   t = [-(fliplr(t)) , 0 , t];
   scalingFactor = (sigmaT * pi^.5)^-.5;
   gaussian = exp(-t.^2 ./ (2 * sigmaT^2));
   sine = exp(1i * 2 * pi * cf .* t);
   kernel = gaussian .* sine;
   kernel = kernel ./ sum(abs(kernel));
   kernel = kernel .* scalingFactor;

   cmw.kernel = kernel'; % transposed to facilitate family table concatenation
   cmw.cf = cf;
   cmw.f1 = f1;
   cmw.f2 = f2;
   cmw.srate = srate;
   cmw.bwo = bwOctaves;
   cmw.bwh = bwHz;
   cmw.sigmaF = sigmaF;
   cmw.ratioCFtoSigmaF = cmw.cf / cmw.sigmaF;
   cmw.sigmaT = sigmaT;
   cmw.nCycles = nCycles;
   cmw.nSigmaT = nSigmaT;
   cmw.length = length(cmw.kernel);
   cmw.halfLength = floor(cmw.length / 2);
   cmw.t = t'; % transposed to facilitate family table concatenation
   cmw.duration = cmw.length / srate;
   cmw.scalingFactor = scalingFactor;
   cmw.mean = mean(kernel);
   cmw.absMean = mean(abs(kernel));
   cmw.lead = kernel(1);
end

function dataConvolved = RTT_Conv(data,cmw)
   % convolves data in time domain with complex morlet wavelet
   [numChans , numPnts] = size(data);
   % pad each end with mirror
   dataPadded = padarray(data,[0,cmw.halfLength],'symmetric','both');
   % pad trail end for integer windowing
   kernelReps = ceil(size(dataPadded,2) / cmw.length);
   dataPadded = padarray(dataPadded,[0,kernelReps * cmw.length - size(dataPadded,2)],'symmetric','post');
   % replicate kernel to same size as dataPad
   kernel = repmat(cmw.kernel',[1,kernelReps]);
   % initialize timeInterval and dataConvolved
   timeInterval = 1 : cmw.length : size(dataPadded,2);
   dataConvolved = complex(NaN(size(data)));
   % perform convolution
   for iStep=0:cmw.length-1
      dataPaddedTemp = circshift(dataPadded,-iStep,2) .* kernel;
      dataPaddedWindowed = reshape(dataPaddedTemp,numChans,cmw.length,kernelReps);
      dataPaddedSum = sum(dataPaddedWindowed,2,'omitnan'); % ignoring NaN
      dataConvolved(:,timeInterval + iStep) = dataPaddedSum;
   end
   % remove padding
   dataConvolved = dataConvolved(:,cmw.halfLength + 1 : cmw.halfLength + numPnts);
end

function KR = RTT_KernelPCA(KU,eegData)
   eegData = double(eegData);
   nVerts = size(KU,1) / 3;
   assert(nVerts == floor(nVerts),'unconstrained kernel not evenly divisible by 3')
   SS = KU * eegData;
   TM = zeros(nVerts, nVerts * 3);
   TMtemp = cell(nVerts,1);
   for i=1:nVerts
      SStemp = SS((i-1)*3+1:i*3,:);
      [eigenVectors, eigenValues] = eig(SStemp * SStemp');
      [~, sortIndices] = sort(diag(eigenValues),'descend');
      TMtemp{i,1} = eigenVectors(:,sortIndices(1))';
   end
   for i=1:nVerts
      TM(i,(i-1)*3+1:i*3) = TMtemp{i};
   end
   KR = TM * KU;
end

function [dataWindowed , numWindows] = RTT_WindowData(data,windowLengthPoints)
   [numChans , numPoints] = size(data);
   numWindows = floor(numPoints / windowLengthPoints);
   numPointsTrunc = numWindows * windowLengthPoints;
   dataWindowed = reshape(data(:,1:numPointsTrunc),numChans,windowLengthPoints,numWindows);
end

function RA = RTT_Orthog(SS)
   % neglible value added to prevent underflow/infinite values from log(0)
   threshOffset = 1e-6;
   % SScda = source space conjugate divided by absolute value, precalculating this saves time in vertex loop
   SScda = conj(SS) ./ abs(SS); 
   % PEX = raw power envelopes
   PEX = SS .* conj(SS);
   % take log to render values more normal
   PEX = log(PEX + threshOffset);
   % precalculate subtracted mean for faster correlation computation
   PEX = PEX - mean(PEX,2);
   % convert to single precision to be gpu practical
%    SS = gpuArray(single(SS));
%    SScda = gpuArray(single(SScda));
%    PEX = gpuArray(single(PEX));
%    RA = gpuArray(single(NaN(size(SS,1))));
   SS    = single(SS);
   SScda = single(SScda);
   PEX   = single(PEX);
   RA    = single(NaN(size(SS,1)));
   
   
   for iVertex=1:size(SS,1)
      % orthogonalize power envelopes with respect to vertex
      PEO = (imag( SS .* SScda(iVertex,:) )).^2;
      % take log to render values more normal
      PEO = log(PEO + threshOffset);
      % correlate, algorithm below faster than corr function
      PEXv = PEX(iVertex,:);
      PEO = PEO - mean(PEO,2);
      RA(:,iVertex) = sum(PEO.*PEXv,2) ./ ((sum(PEO.^2,2)).^.5 .* (sum(PEXv.^2,2)).^.5);
   end
   RA = gather(RA);
end

function [t,p] = RTT_Tstat1(X,muX,tail)
   meanX = mean(X,2,'omitnan');
   stdX = std(X,[],2,'omitnan');
   nX = size(X,2);
   meanDiff = meanX - muX;
   stanError = stdX ./ sqrt(nX);
   df = nX - 1;
   t = meanDiff ./ stanError;
   if strcmpi(tail(1),'B')
      p = 2 * tcdf(-abs(t), df);
   elseif strcmpi(tail(1),'R')
      p = tcdf(-t, df);
   elseif strcmpi(tail(1),'L')
      p = tcdf(t, df);
   end
end
