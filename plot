ft_defaults; 

% Load the Brainnetome (or other) atlas separately 
% Define its filepath, then use ft_read_atlas

% Create variable 'var' with first column as parcel numbers (or other index)
% And second column as the dependent variable values

region_numbers = var(:, 1);  % First column with numbers from 1 to 245
region_values = var(:, 2);  % Second column with values per number

% Step 2: Create a source structure
source = [];
[x, y, z] = ndgrid(1:atlas.dim(1), 1:atlas.dim(2), 1:atlas.dim(3));
coords = [x(:), y(:), z(:)];
coords_homogeneous = [coords, ones(size(coords, 1), 1)]';  % Convert to homogeneous coordinates
coords_transformed = atlas.transform * coords_homogeneous;
source.pos = coords_transformed(1:3, :)';  % Extract x, y, z coordinates
source.dim = atlas.dim;
source.inside = atlas.tissue > 0;  % Mask of inside tissue
source.pow = zeros(prod(atlas.dim), 1);  % Initialize with zeros

% Assign values to the correct regions
for i = 1:length(region_numbers)
    region_index = region_numbers(i);
    value = region_values(i);
    source.pow(atlas.tissue == region_index) = value;
end

% Reshape source.pow to match the volume dimensions
source.pow = reshape(source.pow, atlas.dim);

% Ensure source.inside is correctly defined
source.inside = reshape(source.inside, prod(source.dim), 1);

% Step 3: Plot the data using ft_sourceplot
cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim = 'auto';  % Adjust as needed
cfg.opacitylim = 'auto';   % Adjust as needed
cfg.opacitymap = 'rampup';

% Change the colormap to 'cool'
cfg.colormap = cool;

% Add a title to the plot
cfg.title = 'Brain Regions with Values';

% Plot using 'slice' method
ft_sourceplot(cfg, source);

% Optionally, plot using a different method like 'ortho'
cfg.method = 'ortho';
ft_sourceplot(cfg, source);
