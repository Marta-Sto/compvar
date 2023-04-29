function source_setup_2()
% Setup function for the source power estimation per participant.
% Initialises and defines configurations to later be used in the run
% function that's run on the server.

%% Set up path and add Fieldtrip
% check whether the path is local or on the server and set up path

if ismac % for Mac
    
    fprintf('Hello, we are in Mac');
    
    basepath = '/Volumes/LNDG/Projects/HCP/'; % MEG HCP data
    backend = 'local';
    compile = 'no';
    
    addpath(basepath);
    
    ff = filesep; 
    
    % add Fieldtrip to path
    fieldtrip_dir = '/Users/stojanovic/Documents/Projects/MaxPlanck/Toolboxes/fieldtrip-20220104/'; % on Marta's local Mac
    addpath(fieldtrip_dir);
    % addpath(genpath('/Volumes/LNDG/Projects/HCP/Sourceproject/fieldtrip-20220104/'));
    
    ft_defaults;
   
    % add megconnectome functions                                        
    addpath(genpath('/Volumes/LNDG/Projects/HCP/megconnectome-master/'));

    % basic paths for MEG data
    MODIN   = fullfile(basepath, 'MEG/');
    MODOUT  = fullfile(basepath, 'SOURCEDATA/');
    

else % for TARDIS (remote server)
    
    fprintf('Hello, we are on TARDIS');
    
    % set path
    addpath('/home/mpib/stojanovic/MEG/megconnectome-master/');
    
    basepath = '/home/mpib/stojanovic/';
    addpath(basepath);
    
    % add relevant paths
    addpath('/home/mpib/stojanovic/fieldtrip-20220104/');
    ff = filesep; 
   
    
    ft_defaults;
    
    % basic paths for MEG
    MODIN   = fullfile(basepath, 'MEG/');
    MODOUT  = fullfile(basepath,'SOURCEDATA/');
   
    backend = 'slurm';
    compile = 'no';
    
    overwrite = 1;
    
end

% specify resource allocation on TARDIS
stack = 1;
timreq = 240;      % in minutes per run
memreq = 8*1024^3; % 4GB memory requirement


%check if dir exists, if not, create it
checkdir = exist(MODOUT); % returns 7 if dir is an existing folder
                          % only set up directory if it isn't already one
% if checkdir ~= 7
%     mkdir(MODOUT) % set up dir if 7 is not returned
% end


%% List HCP subjects

% sub = {'100307'; '102816'; '105923'; '106521'; '108323';'109123'; 
%        '111514'; '112920'; '113922'};
% with current script, apparently unable to find several of the subjects

sub = {
    '100307'; '102816'; '105923'; '106521'; '108323'; '109123'; '111514'; 
    '113922'; '116726'; '125525'; '133019'; '140117'; '146129'; '149741'; 
    '153732'; '156334'; '158136'; '162026'; '166438'; '172029'; '175540'; 
    '177746'; '181232'; '185442'; '187547'; '189349'; '191033'; '191437'; 
    '191841'; '192641'; '195041'; '204521'; '205119'; '212318'; '221319'; 
    '233326'; '250427'; '255639'; '287248'; '293748'; '352132'; '352738'; 
    '406836'; '433839'; '512835'; '559053'; '568963'; '599671'; '601127'; 
    '660951'; '662551'; '665254'; '680957'; '715950'; '725751'; '783462'; 
    '814649'; '825048'; '872764'; '877168'; '898176'; '912447'; '917255'}; 

%% Make configurations (cfg) list

% initialise configurations per subject
cfg         = [];
cfglist     = {};

% loop over the configurations across subjects (configurations based on
% what's different between subjects)

for isub = 1:length(sub)
    
    ff = filesep; % file seperator (different for Mac or Windows)
    
    % subject data
    subfolder  = [MODIN sub{isub}];
    restfile   = [subfolder ff [sub{isub} '_' 'MEG_3-Restin_rmegpreproc.mat']];
    anatfile   = [subfolder ff [sub{isub} '_' 'MEG_anatomy_transform.txt']];
    headfile   = [subfolder ff [sub{isub} '_' 'MEG_anatomy_headmodel.mat']]; % head model
    smodfile   = [subfolder ff [sub{isub} '_' 'MEG_anatomy_sourcemodel_3d8mm.mat']]; % source model
    
    % configure input data and output
    cfg.MODIN   = MODIN;  % input
    cfg.MODOUT  = MODOUT; % output
 
    cfg.restfile = restfile; % input file
    cfg.outfile  = [MODOUT sub{isub} 'source_parcel.mat']; % output file
    cfg.anatfile = anatfile;
    cfg.headmod  = headfile;
    cfg.smod     = smodfile;
    cfg.subjno   = sub(isub); % subject number
    
    if ~exist(cfg.outfile, 'file') || overwrite
        cfglist = [cfglist cfg];
    end
end

%% Shuffle configurations and check

% shuffle list of configs to run
cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

% print out the number of data files
fprintf('The source estimation for %d cfgs\n', length(cfglist))

if strcmp(backend, 'slurm')
    options = '-D. -c2'; % --gres=gpu:1
else
    options =  '-l nodes=1:ppn=3';
end

% conduct check
mkdir('~/qsub');


%% Run source estimation and parcellation

if strcmp(compile, 'yes')
    % compile the function to be executed in a batch on the server
    % function called has to live in the same folder as the setup function 
    
    fun2run = qsubcompile({@source_script_run}, 'toolbox', {'signal', 'stats'}); 
else
    fun2run = @source_script_run;                               
end 

% check how this is run
if strcmp(backend, 'local')
    cellfun(fun2run, cfglist)
    return
end

% prepare job list, qsubcellfun applies the source script to each element

qsubcellfun(fun2run, cfglist, 'memreq', memreq, 'timreq', timreq*60, 'stack', stack, ...
    'StopOnError', true, 'backend', backend, 'options', options);

end 

