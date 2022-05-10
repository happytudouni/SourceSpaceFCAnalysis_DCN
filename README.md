These scripts and the EEG Data are part of Xie, W., Toll, R., & Nelson, C.A. (In revision). EEG Functional Connectivity Analysis in the Source Space. Developmental Cognitive Neuroscience. 

# Background
In the current tutorial and manuscript we aimed to demonstrate recently developed pipelines to conduct EEG FC analysis in the source space focusing on phase-phase synchrony (PPS) (Xie et al., 2019a, b, Dev Sci and BMC Med) and and amplitude-amplitude correlation (AAC) (Toll et al., 2020, AMJ; Zhang et al., 2021, Nat Biomed Eng) respectively. The two pipelines were developed by two different research labs -- the [Nelson Lab](https://www.childrenshospital.org/research/labs/nelson-laboratory) at Harvard and the [Etkin Lab](https://www.etkinlab.stanford.edu/) at Stanford. You may also find the information about the first and second authors of the current tutorial at [WX at PKU](https://www.psy.pku.edu.cn/english/people/faculty/professor/wanzexie/index.htm) and [RT at UTSW Med Center](https://profiles.utsouthwestern.edu/profile/185228/russell-toll.html). 
For the sake of brevity, the two pipelines were referred to as the “pl_pps” and “pl_aac” in the following sections. 



# To use pl_pps, please do the followings
1. Download the SourceSpaceFC folder 

    You should see the following files in your folder:

    <img src="https://tva1.sinaimg.cn/large/e6c9d24egy1h22a0trktuj213y08yaai.jpg" alt="Screen Shot 2022-05-09 at 15.01.45" style="zoom:50%;" />

    - Go to the "Data" folder. 
        - Example data have already been uploaded to the "Age12mos" and "Age36mos" folders. You should be able to see them once you open the folders. 
        - To download all the data used in the manuscript, open "Data download.md" file, and then you can download the data following the instructions and have the EEG data for different ages saved in the local "Age12mos" and "Age36mos" folders.
    - Go to the "Sourcemodels" folder. Download the source models and filters following the instructions in the "Download sourcemodels and efilters.md" document and have the files saved in your local "Sourcemodels" folder.
    - Go to the "Programs" folder. Open or run "RunMe_example.m" as an example for how to use the scripts.
        - Make sure to download and install **EEGLAB and Fieldtrip** toolboxes and add them to Matlab search path. We have the programs tested with the following Fieldtrip and EEGLAB versions:
        - [Both fieldtrip-20210409 and fieldtrip-20180415 ](https://www.fieldtriptoolbox.org/download/)
        - [eeglab2021.1](https://eeglab.org/download/)
        - To change the parameters and methods used for source space FC analysis, please see the following scripts as an example. Please type "help SourceSpaceFCanalysis_PPS" in your Matlab command window for more information.

     ```matlab
         help SourceSpaceFCanalysis_PPS
     % SourceSpaceFCanalysis_PPS performs source-space functional connectivity analysis.
     % This program was written based on the example data in Fieldtrip preprocessing format.
     % Use as
     %   SourceSpaceFCanalysis_PPS(cfg, participantnumber)
     % Inputs:
     %   cfg: The configuration (structure) that contains various parameters:
     %   cfg.foi: frequency of interest. The following boundaries are arbitrarily defined based on the literature;
     %            They can be changed to other values. Please read Xie et al.(under revision)_DCN for more information;
     %            =  'theta' (run the analysis for the theta band)
     %                   .12mos:[3 6], which means 3 to 6Hz
     %                   .36mos:[3 7], which means 3 to 7Hz
     %            =  'alpha' (defualt)
     %                   .12mos:[5 10]
     %                   .36mos:[6 11]
     %            =  'beta'
     %                   .[11 22]
     %            =  'beta'
     %                   .[22 45]
     %   cfg.fcmethod: functional connectivity method. Our program calls ft_connectivityanalysis;
     %            Thus, the majority of the fc options in Fieldtrip should work. 
     %            =   'wpli'(default)
     %                   .weighted phase lag index 
     %            =   'imag'
     %                   .imaginary part of the coherency 
     %   cfg.atlastype: brain atlas used for segmentation.
     %            Only the LPBA40 atlas was provided as an example for this tutorial.
     %            Please find more atlases from John E. Richards' website: https://jerlab.sc.edu/projects/neurodevelopmental-mri-database/
     %            =   'LPBA'(default)
     %                   .lpba40 brain atlas 
     %   cfg.replacearg: whether to replace the existing output file.
     %            =   0 (default)
     %                   .Do not replace the existing file. Return if the file already exists
     %            =   1 
     %                   .Replace the existing file
     %   cfg.parcmethod: How to parcellate the voxels into ROIs.
     %            =   'centroid' (default)
     %                   .Use the voxels near the centroid for each ROI
     %            =   'average' 
     %                   .Calculate the average across all voxels for each ROI
     %            =   'PCA' 
     %                   .Do PCA with all voxels in the ROI and use the first component;
     %   cfg.methodtype: Source localization inverse solution.
     %            =   'eloreta' (default)
     %                   .Use eLORETA
     %            =   'mne' 
     %                   .Use minimum norm estimation
     %   cfg.gridresolution: Source localization "spatial resolution".
     %            This is also an arbitrary value; however, you would need to re-generate the source models if you use other values;
     %            Please follow this link to create new models https://www.fieldtriptoolbox.org/tutorial/sourcemodel/
     %            =   '6mm' (defualt)
     %   cfg.plotarg: whether do plots for certain steps while the program is running.
     %            =   0 (default)
     %                   .No plot. Thanks.
     %            =   1 
     %                   .I like plots.
     %   participantnumber: an arbitrary participant number. Used for loading and saving the data.
     ```

2. Run the source-space FC analysis for all participants

     > type ```edit RunAnalysis4All``` in the command window
     > change the parameters for the SourceSpaceFCanalysis_PPS program from line 31 to 40
     > type the name of the program in the command window and hit enter

     ```matlab
     RunAnalysis4All
     ```



# To use pl_aac, please do the followings

1. Download all the files (.mat & .m) in the pl_aac folder.
2. Create a folder with an arbitrary name (e.g., 'pl_aac') and a few subfolders as shown in the following picture:

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gyaqbo972aj30as09a3yk.jpg" alt="Screen Shot 2022-01-12 at 11.27.14 AM" style="zoom:50%;" />

- Save the kernels and .mat files downloaded from 'Sourcemodels/pl_aac_sourcemodels/' to the 	'ADMIN' folder, and download the preprocessed data from "pl_aac_preprocessed" and save them to the 'FINAL' and 'OPEC' folders.

3. The Matlab code should be executed as follows: 
   a. Make sure **Brainstorm** is installed and added to the search path because pl_aac relies on functions in Brainstorm.

   b. on Line 42 in ```babyPublishOrthog.m```, define the 'folderBase' as the directory of the 'pl_aac' folder.

   c. Run the ```babyPublishOrthog.m``` by typing its name ```babyPublishOrthog``` in the command window. 

- Note: The later connectivity steps between pl_aac and pl_pps are not dependent on the earlier source differences.
- Note: We also provide a program "SourceSpaceFCAnalysis_AAC.m" to calculate orthogonalized power correlation (i.e., AAC) with the head models and functions used in pl_pps (e.g., FEM model and Fieldtrip function). The ```OrthogonalPowCorr.m``` was written for this purpose, and it is called by the "SourceSpaceFCAnalysis_AAC" program.  
- Note: The Signal Processing Toolbox (SPT) of Matlab is recommended to be installed to ensure all the programs would run smoothly. The functions in Fieldtrip and Brainstorm might use functions in the SPT, and the customized programs written by us also adopt SPT functions, e.g., the "bandpass.m" used to filter data in "SourceSpaceFCAnalysis_AAC" is a SPT function.  

# Preprocessing programs

In this tutorial we are also sharing the preprocessing programs (mostly automatic).
1. The preprocessing programs used for the pl_pps pipeline were saved in the pl_pps_preproc folder: ```BaselineDataProcessing_pps.m``` and ```DefineArtificialICAs_AdjustAndSASICA.m```. These are the programs we modified from the ones used for Xie et al.(2019), BMC Medicine. The parameters in these programs can be tweaked, and please email WX (wanze.xie@pku.edu.cn) if you have questions about pediatric EEG data preprocessing.
2. The proprocessing program used for the pl_aac pipeline was named as ```babyPublishPreProc```. This program was written based on the preprocessing procedures in Toll et al. (2020), AJP. You may also reach out to RT for questions about this program. 

# Plotting

The brain FC figures depicting the results from using pl_pps were made with [Surf Ice](https://www.nitrc.org/projects/surfice/) software. 
1. Download and install Surf Ice with the link above. 
2. Load a brain image, e.g., "MNI152_2009" by clicking on File --> Open --> Surf Ice --> sample --> mni152_2009.mz3
3. Load a node file, e.g., File --> Open --> BrainNet --> LPBA40.node
4. Load an edge file, e.g., File --> Open --> BrainNet --> LPBA40.edge

The brain FC figures for pl_aac were made with **Brainstorm**

# Intermediate outputs 

In this section, we use one participant's data as an example to show the some of the intermediate outputs from the SourceSpaceFCAnalysis_PPS.m. program.
1. Load the EEG_FT data and run ft_timelockanalysis to change the format of the data from "preprocessing" to "timelock".
```matlab
line 57 to 64
%load the EEG data 
eegfilename = ['Experiment 1 ...']; 
load ([eegfilepath eegfilename],'EEG_FT'); 
% change the data format from preprocessing to timelock format
cfg = [];
cfg.keeptrials = 'yes';
cfg.covariance = 'yes';
EEG_FT = ft_timelockanalysis(cfg,EEG_FT);
```
EEG_FT has the following structure. 

<img src="https://user-images.githubusercontent.com/45924665/147401599-6cb8bfbe-ba07-44a3-9214-c1c3324ecc0c.png" alt="image" style="zoom:50%;" />

2. Cortical source reconstruction --> "sourcedata" that has the following structure
```matlab
from line 153 to line 173
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
```
<img src="https://user-images.githubusercontent.com/45924665/147401633-9cb312c8-c5a4-468c-b90a-6c784dd79576.png" alt="image" style="zoom: 67%;" />

3. Parcellation of the dipoles across the brain into ROIs -->  "roitrialdata" that has the following structure
```
line 175 to line 273
```
<img src="https://user-images.githubusercontent.com/45924665/147401744-6dcbee45-63a9-4f0c-a53e-d690470636e1.png" style="zoom:50%;" />

```rpt_chan_time``` stands for the "number of trials/epoch" x "48 LPBA ROIs" x "time samples (250Hz x 2 seconds)"

4. Frequency analysis of the source-space ROI data --> "EEG_Freq" that has the following structure
```matlab
line 295 to line 303
%% Frequency analysis
    % A few options for cfg.method in Fieldtrip: {'coh', 'csd', 'wpli', 'wpli_debiased', 'plv','imag'}
        cfg            = [];
        cfg.output     = 'fourier';
        cfg.method     = 'mtmfft';
        cfg.foilim     = bpfreq; % 6 to 11 Hz in this example
        cfg.keeptrials = 'yes';
        cfg.taper      = 'hanning';
        EEG_freq       = ft_freqanalysis(cfg, EEG_pm);
```
![](https://user-images.githubusercontent.com/45924665/147401832-ca106856-7b7f-4e21-93ea-b2e2d7a55edc.png)

5. Functional connectivity analysis --> "source_conn" that has the following structure
```matlab
line 305 to line 312
%% FC analysis
        cfg = [];
        F;
        if strcmp(fcmethod,'imag');
            cfg.method = 'coh';
            cfg.complex = 'absimag';
        end
        source_conn = ft_connectivityanalysis(cfg, EEG_freq);
```
![](https://user-images.githubusercontent.com/45924665/147401893-0de73670-3982-4f5c-a480-592af5e26cb3.png) 



# Other programs

- ```textprogressbar.m```   This is a prgram written by Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com) to show a progress bar during cortical source reconstruction. It is not a necessary piece for the source-space FC programs.

- ```GenerateTxtFilesForStats.m```  This program would load the output FC matrix for all participants in one age group (e.g., 12 or 36 mos) and calcuate the average FC between ROIs in the four major lobes, i.e., frontal, temporal, parietal and occipital lobes. 

  - To use this program, do the following:

    ```matlab
    % The following 'cfg' settings are to calculate the average FC in the theta band for the 12 mos group, using the FC outputs from using a combination of the following parameters: 'AAC' for FC calculation, 'average' for ROI parcellation, and 'eLORETA' for source reconstruction.
    cfg = [];
    cfg.foi              = 'theta'; % 'alpaha', 'beta', 'gamma'    
    cfg.fcmethod         = 'AAC';  % or 'wPLI'      
    cfg.replacearg       = 1; % 1: to replace previously generated file; 0: no  replacement         
    cfg.parcmethod       = 'average';  % 'centroid', 'PCA'   
    cfg.methodtype       = 'eloreta';  % 'MNE'         
    cfg.age              = 12;    
    [fcmatrix_all] = GenerateTxtFilesForStats(cfg);
    ```

  - The output of this program is a text file with the participant numbers, age, frequency band, and the FC for different pairs of lobes (e.g., *FF* refers to the average FC between ROIs within the frontal lobe; *FO* refers to the average FC between ROIs in the frontal and those in occipital lobes; and *lpba* is the average FC between all ROIs, i.e., whole-brain FC).

  - This program will also plot the average adjacent matrix across participants.

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gy8mj93nb7j30wy0u0k21.jpg" style="zoom: 50%;" />

