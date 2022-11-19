function source_setup()
% SOURCE SETUP Setup function for the source power estimation per participant.
% Initialises and defines configurations to later be used in the source_run
% function run on the server.

%% Set up paths and add Fieldtrip
% check whether the path is local or on the server

% for Mac
if ismac
    
    % set basepath
    basepath = '/Volumes/LNDG/Projects/HCP/Sourceproject/'; % for MEG HCP data
    backend = 'local';
    compile = 'no';
    
    ff = filesep; 
    
    % add fieldtrip and relevant paths
    fieldtrip_dir = '/Users/stojanovic/Documents/Projects/MaxPlanck/Toolboxes/fieldtrip-20220104/'; % on Marta's local Mac
    addpath(fieldtrip_dir);
    % addpath(genpath('/Volumes/LNDG/Projects/HCP/Sourceproject/fieldtrip-20220104'))                               
    addpath(genpath('/Volumes/LNDG/Projects/HCP/megconnectome-master/')); % add megconnectome functions 

    % basic paths for MEG data
    MODIN   = fullfile(basepath, 'MEG/');
    MODOUT  = fullfile(basepath, 'SOURCEDATA/');
    ATLAS   = fullfile('/Users/stojanovic/Downloads/BNA_MPM_thr25_1.25mm.nii');
    
% for TARDIS
else
    % set basepath
    basepath = '/home/mpib/stojanovic/';
    
    % add relevant paths
    addpath(basepath, 'fieldtrip-20220104/'); % fieldtrip
    addpath(basepath, 'qsub_tardis_slurmpreview/');
    addpath(genpath('/home/mpib/stojanovic/megconnectome-master/')); % megconnectome functions
    % addpath(genpath('/home/mpib/stojanovic/BrainnetomeAtlasViewer-2')); % brainnetome atlas
    
    ft_defaults;
    
    % basic paths for MEG
    MODIN   = fullfile(basepath, 'MEG/');
    MODOUT  = fullfile(basepath,'SOURCEDATA/'); % folder in which to store source data
    ATLAS   = fullfile(basepath, 'BrainnetomeAtlasViewer-2/'); % atlas
    
    backend = 'slurm';
    compile = 'no';
    
    overwrite = 1;
end

% specify resource allocation on TARDIS
stack = 1;
timreq = 240; % in minutes per run
memreq = 8*1024^3; % 4GB memory requirement


%check if dir exists, if not, create it
checkdir = exist(MODOUT); % returns 7 if dir is an existing folder
                          % only set up directory if it isn't already one
if checkdir ~= 7
    mkdir(MODOUT) % set up dir if 7 is not returned
end

%% List HCP subjects - add all once all data is on server
% make sure all data is there before running the function
% sub = {'100307'; '102816'; '104012'; '105923'; '106521'; '108323';'109123'; 
       % '111514'; '112920'; '113922'};

sub = {'100307'};

%% Make configurations (cfg) list
% initialise configurations per subject
cfg         = [];
cfglist     = {};

% loop over the configurations across subjects 
% configurations are based on what's different between subjects

for isub = 1:length(sub)
    
    % configure input data and output
    cfg.MODIN   = MODIN; % input
    cfg.MODOUT  = MODOUT; % output
    cfg.ATLAS   = ATLAS;
    
    ff = filesep; % file seperator (different for Mac or Windows)
    
    % subject data
    subfolder  = [MODIN sub{isub} ff sub{isub} ff 'Restin' ff 'rmegpreproc'];
    restfile   = [subfolder ff [sub{isub} '_' 'MEG_3-Restin_rmegpreproc.mat']];
    anatfolder = [MODIN sub{isub} ff sub{isub} ff 'anatomy'];
    anatfile   = [anatfolder ff [sub{isub} '_' 'MEG_anatomy_transform.txt']];
    headfile   = [anatfolder ff [sub{isub} '_' 'MEG_anatomy_headmodel.mat']]; % head model
    smodfile   = [anatfolder ff [sub{isub} '_' 'MEG_anatomy_sourcemodel_3d8mm.mat']]; % source model
    fmrifile   = [anatfolder ff 'T1w_acpc_dc_restore.nii.gz'];
    % atlasfile  = [ATLAS ff 'Brainnetome_v1.0.2' ff 'Atlas' ff 'atlas.mat'];
    atlasfile  = [ATLAS 'BNA_MPM_thr25_1.25mm.nii'];
        
    cfg.restfile = restfile;  % input file
    cfg.outfile  = [MODOUT sub{isub} 'source_parcel.mat']; % output file
    cfg.anatfile = anatfile;  % anatomy file
    cfg.headmod  = headfile;  % head model
    cfg.smod     = smodfile;  % source model
    cfg.fmri     = fmrifile;  % fmri indir
    cfg.atlas    = atlasfile; % atlas file
    cfg.subjno   = sub(isub); % subject number
    
    if ~exist(cfg.outfile, 'file') || overwrite
        cfglist = [cfglist cfg];
    end
end

%% Shuffle configurations and check
% shuffle list of configs to run
cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

% print out the number of data files
fprintf('The source estimation for %d subject(s)\n', length(cfglist))

if strcmp(backend, 'slurm')
    options = '-D. -c2'; % --gres=gpu:1 % check later these are correct
else
    options =  '-l nodes=1:ppn=3';
end

% conduct check
mkdir('~/qsub'); cd('~/qsub');

%% Run source estimation on parcellated data
if strcmp(compile, 'yes')
    % compile the function to be executed in a batch on the server
    % function called has to live in the same folder as the setup function 
    % for the handle to work
    
    fun2run = qsubcompile({@source_run}, 'toolbox', {'signal', 'stats'}); 
else
    fun2run = @source_run;    % fun2run takes all prior inputs and returns 
                              % the parcellated source estimation per subject                                   
end 

% check how this is run
if strcmp(backend, 'local')
    cellfun(fun2run, cfglist)
    return
end

% prepare job list, where qsubcellfun applies the source script to each 
% element of a cell array

qsubcellfun(fun2run, cfglist, 'memreq', memreq, 'timreq', timreq*60, 'stack', stack, ...
    'StopOnError', true, 'backend', backend, 'options', options);
