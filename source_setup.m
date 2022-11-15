function source_setup()
% Setup function for the source power estimation per participant.
% Initialises and defines configurations to later be used in another
% function run on the server.

%% Set up path and add Fieldtrip
% check whether the path is local (Mac) or on the server and set up path

% for Mac
if ismac
    % using Mountain Duck to connect to the appropriate location on the server
    % basepath = '/Users/stojanovic/Library/Group_Containers/G69SCX94XU.duck/Library/Application_Support/duck/Volumes/tardis.mpib-berlin.mpg.de - SFTP/ScriptFolders';
    
    %basepath = '/Volumes/LNDG/Projects/HCP/'; % for MEG HCP data
    basepath = '/home/mpib/stojanovic/ScriptFolders/';
    backend = 'local';
    compile = 'no';
    
    % add Fieldtrip to path
    fieldtrip_dir = '/Users/stojanovic/Documents/Projects/MaxPlanck/Toolboxes/fieldtrip-20220104'; % on Marta's local Mac
    addpath(fieldtrip_dir);
    
    % run ft_defaults
    ft_defaults
   
    % add megconnectome functions                                        
    addpath(genpath('/Volumes/LNDG/Projects/HCP/megconnectome-master/'));

    % basic paths for MEG data
    MODIN   = fullfile(basepath, 'MEG/');
    MODOUT  = fullfile(basepath, 'SOURCEDATA');
    

% for TARDIS (remote server)
else
    
    % when running from Matlab on TARDIS
    % basepath = 'ScriptFolders/';
    
    % set path and add fieldtrip
    basepath = '/home/mpib/stojanovic/ScriptFolders/';
    addpath([basepath, 'fieldtrip-20220104'])
    ft_defaults;
    
    % basic paths for MEG
    MODIN   = fullfile(basepath,'MEG/');
    MODOUT  = fullfile(basepath,'SOURCEDATA/'); % folder in which to store source data
   
    
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

%% List HCP subjects - add all once all data is loaded
% make sure all data is there before running script
sub = {'100307'; '102816'; '104012'; '105923'; '106521'; '108323';'109123'; 
       '111514'; '112920'; '113922'};

%% Make configurations (cfg) list
% initialise empty list and cell for configurations
cfg         = [];
cfglist     = {};

% loop over the configurations across subjects 
% configurations are based on what's different between subjects

for isub = 1:length(sub)
    
    ff = filesep; % file seperator (different for Mac or Windows)
    
    % subject data
    subfolder  = [MODIN sub{isub} ff 'Restin' ff 'rmegpreproc'];
    restfile   = [subfolder ff [sub{isub} '_' 'MEG_3-Restin_rmegpreproc.mat']];
    anatfolder = [MODIN sub{isub} ff 'anatomy']; % anatomy
    anatfile   = [anatfolder ff [sub{isub} '_' 'MEG_anatomy_transform.txt']];
    transform  = [hcp_read_ascii(anatfile)]; % transformation matrix per participant
    headfile   = [anatfolder ff [sub{isub} '_' 'MEG_anatomy_headmodel.mat']]; % head model
    smodfile   = [anatfolder ff [sub{isub} '_' 'MEG_anatomy_sourcemodel_3d8mm.mat']]; % source model
    % source     = [MODIN sub{isub} ff 'Source' ff];
    
    % configure input data and output
    cfg.MODIN   = MODIN; % input
    cfg.MODOUT  = MODOUT; % output
    
    % configure source data
    % cfg.SOURCE  = source; 
 
    cfg.infile   = restfile; % input file
    cfg.outfile  = [MODOUT sub{isub} 'source_parcel.mat']; % output file
    cfg.anatfile = anatfile;
    cfg.trans    = transform;
    cfg.headmod  = headfile;
    cfg.smod     = smodfile;
    cfg.subjno = sub(isub); % subject number
    
    if ~exist(cfg.outfile, 'file') || overwrite
        cfglist = [cfglist cfg];
    end
end

%% Shuffle configurations and check
% shuffle list of configs to run
cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

% print out the number of data files
% fprintf('The source estimation for %d cfgs\n', length(cfglist))

% if strcmp(backend, 'slurm')
    % options = '-D. -c2'; % --gres=gpu:1 % check later these are correct
% else
    % options =  '-l nodes=1:ppn=3';
% end

% conduct check
% mkdir('~/qsub'); cd('~/qsub');

%% Run source estimation and parcellation
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

% prepare job list, where qsubcellfun applies the source script 
% to each element of a cell array

qsubcellfun(fun2run, cfglist, 'memreq', memreq, 'timreq', timreq*60, 'stack', stack, ...
    'StopOnError', true, 'backend', backend, 'options', options);

