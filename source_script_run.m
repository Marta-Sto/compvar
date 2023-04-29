function source_script_run(cfg)
% Function to run source analyses on data.

%% add relevant paths and define information

datdir = '/home/mpib/stojanovic/SOURCEDATA/';
addpath(datdir);

%% create modout directory
MODOUT = '/home/mpib/stojanovic/SOURCEDATA/';

%check if dir exists, if not, create it
checkdir = exist(MODOUT); % returns 7 if dir is an existing folder
          
if checkdir ~= 7
    mkdir(MODOUT) % set up dir if 7 is not returned
end

%% Pass cfg list
% pass the configurations set up in the setup function

restfile   = cfg.restfile; 
anatfile   = cfg.anatfile;
headfile   = cfg.headmod;
smodfile   = cfg.smod; 
subjno     = cfg.subjno;

%% load relevant data per subject
load(restfile);

% transformation matrix per subject
tmat = hcp_read_ascii(anatfile); 
T = tmat.transform.bti2spm;

% head model 
tmp = load(headfile);
hm = tmp.headmodel;
hm = ft_convert_units(hm, 'mm');

% source model
tmp = load(smodfile);
sm = tmp.sourcemodel3d;
sm = ft_convert_units(sm, 'mm');
disp('Headmodel and sourcemodel ready');

%% preprocess trials
clear trls;
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = 'all';
trls = ft_preprocessing(cfg,data);
disp('Done preprocessing and baseline correcting trials');

%% get covariance of trials
clear tlck_first;
cfg = [];
cfg.covariance  = 'yes';
tlck_first = ft_timelockanalysis(cfg,trls);
disp('Done timelocking trials');

%% get covariance across all 10-second trials
clear tlck_all;
cfg = [];
cfg.keeptrials   = 'yes';
cfg.vartrllength = 0;
cfg.covariance   = 'yes';
tlck_all         = ft_timelockanalysis(cfg,trls);
disp('Done timelocking all trials');

%% also transform headmodel
clear hm_trans;
hm_trans = ft_transform_geometry(T,hm);
disp('Done transforming head model');

%% get reg kappa and run beamformer (LCMV)
clear kappa;
[u,s,v] = svd(((tlck_first.cov)));
d       = -diff(log10(diag(s)));
d       = d./std(d);
kappa   = find(d>5,1,'first');

%% generate leadfield
clear lf_meg;
clear lf_meg_trans;
cfg             = [];
cfg.channel     = tlck_first.label;
cfg.grad        = tlck_first.grad;
cfg.normalizeparam = 1;
cfg.sourcemodel = sm;
cfg.headmodel   = hm;
cfg.reducerank  = 2;
cfg.normalize   = 'yes';
lf_meg       = ft_prepare_leadfield(cfg);
lf_meg_trans = ft_transform_geometry(T,lf_meg);
disp('Done creating leadfield');

%% run source analysis on longer trials
clear source_avg;
cfg                 = [];
cfg.lcmv.lambda     = '5%';
cfg.lcmv.kappa      = kappa;
cfg.method          = 'lcmv';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.weightnorm = 'unitnoisegain';
cfg.headmodel   = hm;
cfg.sourcemodel = lf_meg; % can also visualise using the transformed version to check
source_avg      = ft_sourceanalysis(cfg,tlck_all);
disp('Done averaging source data');

%% read atlas
clear brainnetome;
brainnetome = ft_read_atlas([fieldtrip_dir '/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii']);
disp('Atlas read');

%% interpolate source average onto atlas
clear atlas_lowres;
cfg = [];
cfg.parameter = 'all';
cfg.interpmethod = 'nearest';
atlas_lowres = ft_sourceinterpolate(cfg,brainnetome,source_avg);

% align positions of the atlas to the source average
atlas_lowres.pos = source_avg.pos;
disp('Done interpolating source data onto the atlas');

%% parcellate single trials
clear parc_trls;
cfg = [];
cfg.method  = 'svd'; % can be either svd or pca
cfg.parcellation = 'tissue';
parc_trls = ft_virtualchannel(cfg,tlck_all,source_avg,atlas_lowres);
disp('Done parcellating single trial data');

%% save the parcellated single trial data
% create directory to save the source data (if it doesn't already exist)
if ~exist([datdir 'Source_parc_data' ff]) % if the directory does not exist
    mkdir([datdir 'Source_parc_data' ff])
end

% save the data to the directory
save([datdir 'Source_parc_data' ff 'CompVar_source_parc_trls' '_' cell2mat(subjno)], 'parc_trls',...
        '-V7.3')

end
