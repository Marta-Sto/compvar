% basic MEG source projection rationale / example
% This is certainly not the most advanced thing one could do and we
% probably can step it up a bit . Nevertheless, this is an OK place to
% start.
clc
clear
% add fieldtrip
fieldtrip_dir = '/Users/waschke/Documents/Matlabtoolboxes/fieldtrip-20220104';
addpath(fieldtrip_dir)
ft_defaults
% also add megconnectome functions
addpath(genpath('/Volumes/LNDG/Projects/HCP/megconnectome-master/'));

rootdir = '/Volumes/FB-LIP/LNDG/Projects/HCP/Sourceproject/';
% get filesep depending on OS
ff = filesep;
megdir = [rootdir 'MEG' ff];
subs_avail = dir(megdir);
% as so often, dir will  also report on some not really present files.
% -->Always important to check subs_avail before looping through names.
% Example participant:
exmpl_sub = subs_avail(5).name;
datdir = [megdir exmpl_sub ff];
% load some resting state MEG data, eyes open (4 & 5 are the following scans)
fname = dir([datdir ff 'Rest*' ff 'rmeg*' ff exmpl_sub '*MEG_3' '*.mat']);
load([fname.folder ff fname.name])
% get transformation matrx for given subject
hcp_read_ascii([datdir 'anatomy' ff exmpl_sub '_MEG_anatomy_transform.txt'])
T = transform.bti2spm;


% load headmodel, convert to milimiters (maybe useful later)
headmod = dir([datdir 'an*' ff exmpl_sub '*headmodel.mat']);
tmp = load([headmod.folder  ff headmod.name]);
hm = tmp.headmodel;
hm = ft_convert_units(hm, 'mm');
hm_trans = ft_transform_geometry(T, hm);
hm_trans.coordsys = 'mni';
% note that this is not exactly like in the tutorial. I tried to use
% subject specific head models here but you could try and use the standard
% one like suggested in the tutorial.

bntm = ft_read_atlas([fieldtrip_dir '/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii']);

% load source model
smod = dir([datdir 'an*' ff exmpl_sub '*source*' '3d8*' '.mat']);
tmp = load([smod.folder  ff smod.name]);
sm = tmp.sourcemodel3d;
sm = ft_convert_units(sm, 'mm');
sm_trans = ft_transform_geometry(T, sm);
sm_trans.coordsys = 'mni';

%% apply atlas
cfg = [];
cfg.atlas      = bntm;
cfg.roi        = bntm.tissuelabel;  % here you can also specify a single label, i.e. single ROI
mask           = ft_volumelookup(cfg, sm_trans);

% determine inside positions
sm_trans.inside = false(sm_trans.dim);
sm_trans.inside(mask==1) = true;

%% plot
figure;
ft_plot_mesh(sm_trans.pos(sm_trans.inside,:));

%% load mri
mri = ft_read_mri([datdir 'anatomy' ff 'T1w_acpc_dc_restore.nii.gz']);
mri.coordsys = 'mni';

cfg                = [];
cfg.warpmni   = 'yes';
cfg.template  = sm_trans;
cfg.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri;
sourcemodel_to_use        = ft_prepare_sourcemodel(cfg);

%% filter
% note that we probably would not do this for task data and certianly would
% use a different filter for any event-related design. Here, this might be
% fine.
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 80];
data = ft_preprocessing(cfg, data);
%% get covariance
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = 'all';
cfg.covariance = 'yes';
tlck = ft_timelockanalysis(cfg, data);

%% tlck of all trials
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = 'all';
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
tlck_all = ft_timelockanalysis(cfg, data);

%% generate leadfield
cfg             = [];
cfg.channel     = tlck.label;
cfg.grad        = tlck.grad;
cfg.normalizeparam = 1;
cfg.sourcemodel = sourcemodel_to_use;
cfg.headmodel   = hm_trans;
cfg.reducerank = 2;
cfg.normalize = 'yes';
lf_meg   = ft_prepare_leadfield(cfg);
lf_meg_trans = ft_transform_geometry(T, lf_meg);

%% get reg kappa and run beamformer (LCMV)
[u,s,v] = svd(((tlck.cov)));
d       = -diff(log10(diag(s)));
d       = d./std(d);
kappa   = find(d>5,1,'first');


cfg                 = [];
cfg.lcmv.lambda = '5%';
cfg.lcmv.kappa = kappa;
cfg.method          = 'lcmv';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.weightnorm = 'unitnoisegain';
cfg.headmodel   = hm_trans;
cfg.sourcemodel = lf_meg;
%cfg.rawtrial = 'yes';
source_avg          = ft_sourceanalysis(cfg, tlck);

%% make atlas low res
cfg = [];
cfg.parameter = 'all';
cfg.interpmethod = 'nearest';
atlas_lowres = ft_sourceinterpolate(cfg, bntm, source_avg);
atlas_lowres.pos = source_avg.pos;

%% parcellate single trials
cfg = [];
cfg.method = 'pca';
cfg.parcellation = 'tissue';
parc_trls = ft_virtualchannel(cfg, tlck_all, source_avg, atlas_lowres);
