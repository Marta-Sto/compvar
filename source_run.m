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

restfile   = cfg.restfile; 
anatfile   = cfg.anatfile;
headfile   = cfg.headmod;
smodfile   = cfg.smod;
fmrifile   = cfg.fmri;
atlasfile  = cfg.atlas;   
subjno     = cfg.subjno;
inflated   = cfg.inflated;

% print out current subject number
SUBJ = sprintf(cell2mat(subjno));

disp(cfg.MODIN);
cd(cfg.MODOUT);

ff = filesep;

%% Load data, headmodel and sourcemodel 

% data per subject 
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

%% Filter
% this would not be done for task data or an event-related design. 

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 80];
data = ft_preprocessing(cfg, data);
disp('Done with filtering');

%% Get covariance

cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = 'all';
cfg.covariance = 'yes';
tlck = ft_timelockanalysis(cfg, data);
disp('Done with Timelock');

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
disp('Leadfield ready');

%% get reg kappa and run beamformer (LCMV)
% Source analysis

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
disp('Sourcemodel ready')

%% Combine right and left hemispheres

% combine R and L hemispheres of template
if exist('template_surface_inflated.mat')
else
fname = '/home/mpib/stojanovic/megconnectome-master/template/Conte69.L.midthickness.4k_fs_LR.surf.gii';
% read cortical sheet
inflated = ft_read_headshape({fname strrep(fname, '.L.', '.R.')});
save('template_surface_inflated', cell2mat(subjno), '.mat', 'inflated')
end

%% Project single trials

disp('Preparing to project single trials')
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

% create directory to save the source trial data
if ~exist('/home/mpib/stojanovic/SOURCEDATA/')
    mkdir('/home/mpib/stojanovic/SOURCEDATA/')
end

% -V7.3 takes a very long time, but is required for bigger data structures.
% can use alternatives - some listed in the TARDIS documentation.
save(['/home/mpib/stojanovic/SOURCEDATA' ff 'CompVar_source_trials_folcmv_sub_' cell2mat(subjno)], 'SourceTrialData',...
        '-V7.3')

%% Interpolate with mri

disp('Interpolating')
% read mri file
ind_mri = ft_read_mri(fmrifile);

% read template mri
template_mri = ft_read_mri('/home/mpib/stojanovic/megconnectome-master/template/mni_icbm152_t1_tal_nlin_sym_09a.nii');

mri = ft_volumereslice(cfg,ind_mri);

cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'pow';
source_mri = ft_sourceinterpolate(cfg,source,ind_mri);
sc = ft_volumenormalise(cfg,source_mri);

%% Parcellate

disp('Parcellating')
% % read atlas
atlas = ft_read_atlas(atlasfile);

% interpolate the source data with the atlas
clear cfg
cfg = [];
cfg.interpmethod = 'nearest';
atlas.coordsys = 'mni'; 
cfg.parameter    = 'pow';
% source_parc = ft_sourceinterpolate(cfg, atlas, sc);
source_parc = ft_sourceinterpolate(cfg,sc,atlas); 

% atlas.pos = source_parc.pos;
cfg.parcellation   = 'parcellation';
cfg.template       = 'atlas';
cfg.parameter      = 'all';
cfg.method         = 'mean';
parc_source = ft_sourceparcellate(cfg, source_parc, atlas);

disp('Done parcellating, and almost all done, well done!')

%% Save source parcellated data

disp('Saving parcellated source data')
save(['/home/mpib/stojanovic/SOURCEDATA' ff 'parc_source' cell2mat(subjno)], 'parc_source',...
        '-V7.3')

end
