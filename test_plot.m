%% Script to plot ICC and other values on surface brain plot, where values
%% plotted are stored in the variable 'parcel_values_test'.

% Load Fieldtrip
ft_defaults;

% Load the Brainnetome atlas (assuming you have the atlas file in NIfTI format)
atlas = ft_read_atlas('BNA_MPM_thr25_1.25mm.nii');

% Read the MRI template (this is often needed for interpolation)
mri = ft_read_mri('standard_mri.mat');

% Prepare the atlas to source structure
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
cfg.mri = mri; % Use the standard MRI as reference
source = ft_sourceinterpolate(cfg, atlas, mri);

% Assuming 'parcel_values_test' is a vector with values for each parcel
% Ensure 'parcel_values_test' has the same number of elements as unique parcel labels in the atlas
unique_labels = unique(source.tissue(:)); % Get unique parcel labels

% Check if the length of parcel_values_test matches the number of parcels
if length(unique_labels) ~= length(parcel_values_test)
    error('The length of parcel_values_test does not match the number of unique parcels in the atlas.');
end

% Initialize the source structure for plotting
source.pow = zeros(size(source.tissue));

% Assign values from 'parcel_values_test' to each parcel
for i = 1:length(unique_labels)
    label = unique_labels(i);
    if label == 0
        continue; % Skip if label is 0 (background)
    end
    idx = find(source.tissue == label);
    source.pow(idx) = parcel_values_test(i); % Assign your value to the parcel
end

% Configure the plotting
cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.funcolormap  = 'jet';
cfg.projmethod   = 'nearest';
cfg.surffile     = 'surface_white_both.mat'; % Use a suitable surface file for visualization
cfg.colorbar     = 'yes'; % Add a color bar

% Plot the data
ft_sourceplot(cfg, source);

