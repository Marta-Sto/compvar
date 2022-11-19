function source_run(cfg)
%SOURCE RUN Takes preprocessed HCP MEG data and returns the source
%estimates and source projected data.

% Parcellation obtained based on the Brainnetome atlas.

% To generate activity and power estimates for parcels based on the atlas, 
% the following steps are conducted:
% 1) read the preprocessed HCP MEG data, 
% 2) read the headmodel and the source model, 
% 3) filter
% 4) generate a leadfield, 
% 5) project single trials,
% 6) interpolate and normalize the data, 
% 7) visualise the source projected data, and finally
% 8) obtain estimates of activity per parcel (and save this data).

%% Pass cfg list
% pass the configurations set up in the setup function

% MODIN = cfg.MODIN;
% MODOUT = cfg.MODOUT;

restfile   = cfg.restfile; 
outfile    = cfg.outfile;
anatfile   = cfg.anatfile;
headfile   = cfg.headmod;
smodfile   = cfg.smod;
fmrifile   = cfg.fmri;
atlasfile  = cfg.atlas;   
subjno     = cfg.subjno;

% print out current subject number
SUBJ = sprintf(cell2mat(subjno));

disp(cfg.MODIN);
cd(cfg.MODOUT);

ff = filesep;

%% Load data, headmodel and sourcemodel 

% data 
load([restfile]); % loads data per subject

% transformation matrix
% hcp_read_ascii([anatfile]); % transformation matrix per participant
tmat = hcp_read_ascii([anatfile]);
T = tmat.transform.bti2spm;

% head model 
tmp = load(headfile);
hm = tmp.headmodel;
hm = ft_convert_units(hm, 'mm');

% source model
tmp = load(smodfile);
sm = tmp.sourcemodel3d;
sm = ft_convert_units(sm, 'mm');

%% Filter
% this would not be done for task data or an event-related design. 

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 80];
data = ft_preprocessing(cfg, data);

%% Get covariance

cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = 'all';
cfg.covariance = 'yes';
tlck = ft_timelockanalysis(cfg, data);

%% Generate leadfield

cfg             = [];
cfg.channel     = tlck.label;
cfg.grad        = tlck.grad;
cfg.normalizeparam = 1;
cfg.sourcemodel = sm;
cfg.headmodel   = hm;
cfg.reducerank  = 2;
cfg.normalize   = 'yes';
lf_meg          = ft_prepare_leadfield(cfg);

%% get reg kappa and run beamformer (LCMV)
[u,s,v] = svd(tlck.cov);
d       = -diff(log10(diag(s)));
d       = d./std(d);
kappa   = find(d>5,1,'first');

cfg                 = [];
cfg.lcmv.lambda = '5%';
cfg.lcmv.kappa = kappa;
cfg.method          = 'lcmv';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'no';
cfg.lcmv.weightnorm = 'unitnoisegain';
cfg.headmodel   = hm;
cfg.sourcemodel = lf_meg;
source          = ft_sourceanalysis(cfg, tlck);

%% surface

% combine R and L hemispheres of template
if exist('template_surface_inflated.mat')
    
else
fname = '/home/mpib/stojanovic/megconnectome-master/template/Conte69.L.midthickness.4k_fs_LR.surf.gii';

inflated = ft_read_headshape({fname strrep(fname, '.L.', '.R.')});

save('template_surface_inflated.mat', subjno, 'inflated')

end

%% Project single trials

% main output that we want to save
clear SourceTrialData
SourceTrialData.trial = cell(1, length(data.trial));
SourceTrialData.label = cell(1,1);
SourceTrialData.trialinfo = data.trialinfo;
SourceTrialData.fsample = data.fsample;
SourceTrialData.time = data.time;
NTrials = length(data.trial);

% get trafo-matrix
InVox       = source.inside;
TrafoMatrix = permute(squeeze(cat(3,source.avg.filter{InVox})),[3 2 1]);

% get labels right
SourceTrialData.label = [];
for vx=1:length(find(InVox==1))
    SourceTrialData.label{vx,1} = strcat('VirtCH', num2str(vx));
end

% source-project single trials
for tr = 1:NTrials
    display(['Projecting trial ' num2str(tr)])
    for mor = 1:3
    SourceTrialData.trial{tr}(:,:,mor) = TrafoMatrix(:,:,mor)*data.trial{tr};
    %SourceTrialData.trial{tr}(:,:) = TrafoMatrix(:,:).*data.trial{tr};
    end
end

%% Save source trial data 

% change all save directories
% make directory to save the source projected data
if ~exist('/home/mpib/stojanovic/SOURCEDATA/')
    mkdir('/home/mpib/stojanovic/SOURCEDATA/')
end
% -V7.3 takes a very long time, but is required for bigger data structures.
% can use alternatives - some listed in the TARDIS documentation.
save(['/home/mpib/stojanovic/SOURCEDATA/' ff 'CompVar_source_trials_folcmv_sub_' subjno], 'SourceTrialData',...
        '-V7.3')

% filename = fullfile(datdir, 'Source',  sprintf('%s_source_lcmv', exmpl_sub));
% save(filename, 'source', 'tlckw');

% also save / source as dummy source
save(['/home/mpib/stojanovic/SOURCEDATA/CompVar_dummy_source_folcmv_sub_100307/'], 'source'); 

%% Interpolate and normalise

% read mri for individual subject
ind_mri = ft_read_mri([fmrifile]);

% read template mri
template_mri = ft_read_mri('/home/mpib/stojanovic/megconnectome-master/template/mni_icbm152_t1_tal_nlin_sym_09a.nii');
clear dummy_source

% activation   = alpha_dat.powspctrm;
dummy_source = source;
dummy_source.avg.pow
% dummy_source.avg.pow(dummy_source.inside) = activation;
dummy_source.time = 0;

% before we interpolate source-projected data with the template,
% we need to transform the coordinates.
dummy_source = ft_transform_geometry(T, dummy_source);

% interpolate sources
cfg = [];
cfg.parameter        = 'all';
source_intp.coordsys = 'mni';
cfg.interpmethod     = 'linear';
source_intp = ft_sourceinterpolate(cfg, dummy_source, ind_mri);

% normalize
cfg = [];
cfg.template = '/home/mpib/stojanovic/megconnectome-master/template/mni_icbm152_t1_tal_nlin_sym_09a.nii';
cfg.templatecoordsys = 'mni';

source_intp = ft_volumenormalise(cfg, source_intp);

%% Plot

n_lgth = size(source_intp.inside,1)*size(source_intp.inside,2)*...
    size(source_intp.inside,3);
source_f_surf = source_intp;
source_f_surf.inside = reshape(source_intp.inside,[n_lgth,1]);
source_f_surf.pow = reshape(source_intp.pow,[n_lgth,1]);
source_f_surf.eta = reshape(source_intp.eta,[n_lgth,1]);
%source_f_surf.pos = reshape(source_intp.pos,[n_lgth,1]);
source_f_surf = rmfield(source_f_surf, 'mni');

cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'pow';

cfg.projmethod = 'project';
cfg.renderer   ='opengl';
cfg.sphereradius   = 1;
cfg.maskparameter  = cfg.funparameter;
cfg.surfdownsample = 10;
%cfg.funcolorlim   = [0 max(abs(source_intp.pow))];
cfg.camlight ='yes';
cfg.surffile = 'template_surface_inflated.mat';
cfg.interactive = 'yes';

% view from left
ft_sourceplot(cfg,source_f_surf);

%% Parcellate

% read atlas
% have to read the atlas another way
% brainnetome = ft_read_atlas('template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');
brainnetome = ft_read_atlas(atlasfile);

cfg = [];
cfg.template     = 'brainnetome';
cfg.interpmethod = 'nearest';
cfg.parameter    = 'pow';
cfg.visible      = 'on';
% interpolate the source data with the atlas.
% this is done to align the source data with the atlas regions.
% run on server given it's too large to run locally.

intp = ft_sourceinterpolate(cfg, source_intp, brainnetome);

cfg = [];
cfg.method = 'vertex';
cfg.funparameter = 'pow';
cfg.anaparameter = 'tissue';
cfg.atlas = 'brainnetome';

% plot source based on the parcellated data
ft_sourceplot(cfg, source, brainnetome);

% source parcellate data based on atlas
parc_source = ft_sourceparcellate(cfg, intp, brainnetome);

%% Save source parcellated data



end
