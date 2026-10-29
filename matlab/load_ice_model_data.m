function load_ice_model_data()
% LOAD_ICE_MODEL_DATA Load separately obtained ICE data into the base workspace.
% Set ICE_MODEL_DIR to an authorised model-data folder, or place the six
% required files in matlab/private_ice/ (excluded from version control).

model_dir = getenv('ICE_MODEL_DIR');
if isempty(model_dir)
    model_dir = fullfile(fileparts(mfilename('fullpath')), 'private_ice');
end

% Require an absolute directory so run() cannot change how later paths resolve.
[found, attributes] = fileattrib(model_dir);
if ~found || ~attributes.directory
    error('ICE:MissingModelData', ...
        ['ICE model data directory is missing. Obtain the model separately ' ...
         'and set ICE_MODEL_DIR to its folder, or use matlab/private_ice/. ' ...
         'See EXTERNAL_MODEL.md. Requested folder: %s'], model_dir);
end
model_dir = attributes.Name;

scripts = { ...
    'ukicematlab.m', ...
    'constants.m', ...
    'enginedata.m', ...
    'airdensity.m', ...
    'airtemperature.m', ...
    'aerolookupdataall.m' ...
};

missing = scripts(~cellfun(@(name) isfile(fullfile(model_dir, name)), scripts));
if ~isempty(missing)
    error('ICE:MissingModelData', ...
        'ICE model data is incomplete in %s. Missing: %s. See EXTERNAL_MODEL.md.', ...
        model_dir, strjoin(missing, ', '));
end

% The supplied scripts call one another by name. Check all files before
% changing the path or executing any script, and restore the caller's path.
previous_path = path;
restore_path = onCleanup(@() path(previous_path)); %#ok<NASGU>
addpath(model_dir, '-begin');
for k = 1:numel(scripts)
    script_path = strrep(fullfile(model_dir, scripts{k}), '''', '''''');
    evalin('base', sprintf('run(''%s'')', script_path));
end
end
