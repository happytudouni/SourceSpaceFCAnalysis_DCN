% Run analysis for all subjects;
ages = [12 36];
participantnumbers  = [];
for i = 1:length(ages)
    age = ages(i);
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
end 


for i = 1:length(ages)
    age = ages(i);
    subjects = participantnumbers.(['age' num2str(age)]);
    for ii = 1:length(subjects)
    participantnumber = subjects(ii);
    % functional connectivity analysis  
    fois = {'theta','alpha','beta','gamma'};
    for iii = 1:length(fois)
        cfg            = [];
        cfg.foi        = fois{iii};
        cfg.fcmethod   = 'wpli';    %'wpli','imag','coh' or other methods in ft_fcanalysis
        cfg.parcmethod = 'average'; %'centroid' or 'average'
        cfg.atlastype  = 'LPBA'; 
        cfg.methodtype = 'eloreta'; %'eloreta'  or 'MNE'
        cfg.age        = age;
        cfg.replacearg = 0;
        cfg.plotarg    = 0;
        SourceSpaceFCanalysis_PPS(cfg, participantnumber);
    end
    end % for ii = 1:length(subjects)
end % for i = 1:length(ages)


%% generate txt files
%{
cfg = [];
cfg.age        = 12;
cfg.foi        = 'theta';
cfg.fcmethod   = 'wpli';
cfg.parcmethod = 'centroid';
cfg.methodtype = 'eloreta';
[fcmatrix_all] = GenerateTxtFilesForStats(cfg);

cfg.foi        = 'alpha';
[fcmatrix_all] = GenerateTxtFilesForStats(cfg);
cfg.foi        = 'beta';
[fcmatrix_all] = GenerateTxtFilesForStats(cfg);
cfg.foi        = 'gamma';
[fcmatrix_all] = GenerateTxtFilesForStats(cfg);


cfg.age        = 36;
cfg.foi        = 'theta';
[fcmatrix_all] = GenerateTxtFilesForStats(cfg);
cfg.foi        = 'alpha';
[fcmatrix_all] = GenerateTxtFilesForStats(cfg);
cfg.foi        = 'beta';
[fcmatrix_all] = GenerateTxtFilesForStats(cfg);
cfg.foi        = 'gamma';
[fcmatrix_all] = GenerateTxtFilesForStats(cfg);
%}
%{
for i = 1:length(ages)
    age = ages(i);
    subjects = participantnumbers.(['age' num2str(age)]);
    for ii = 1:length(subjects)
    participantnumber = subjects(ii);
    % functional connectivity analysis  
    fois = {'theta','alpha','beta','gamma'};
    for iii = 1:length(fois)
        cfg            = [];
        cfg.foi        = fois{iii};
        cfg.parcmethod = 'average'; %'centroid' or 'average'
        cfg.atlastype  = 'LPBA'; 
        cfg.methodtype = 'MNE';     %'eloreta' or 'MNE'
        cfg.age        = age;
        cfg.replacearg = 0;
        cfg.plotarg    = 0;
        SourceSpaceFCanalysis_AAC(cfg, participantnumber);
    end
    end % for ii = 1:length(subjects)
end % for i = 1:length(ages)
%}