function convert_for_python(inputFile)
% CONVERT_FOR_PYTHON Converts Simulink logsout data to a Python-friendly struct.
%
%   Usage:
%       convert_for_python()          % Prompts for file
%       convert_for_python(filePath)  % Uses provided file
%
%   Output:
%       Saves a .mat file with suffix "_py.mat" in the same directory.
%       The .mat file contains a struct 'data' where each field corresponds
%       to a signal from logsout.

% --- 1. Handling Input ---
if nargin < 1 || isempty(inputFile)
    [fileName, folderPath] = uigetfile('*.mat', 'Select Simulation Result File');
    if fileName == 0; return; end
    inputFile = fullfile(folderPath, fileName);
end

fprintf('Processing: %s\n', inputFile);
loadedData = load(inputFile);

% --- 2. Locate logsout ---
if isfield(loadedData, 'out') && isprop(loadedData.out, 'logsout')
    logs = loadedData.out.logsout;
elseif isfield(loadedData, 'logsout')
    logs = loadedData.logsout;
else
    error('Could not find "out.logsout" or "logsout" in the file.');
end

% --- 3. Extract and Flatten Data ---
data = struct();
numElements = logs.numElements;

for i = 1:numElements
    elem = logs.get(i);
    name = elem.Name;

    fprintf('  Processing signal %d/%d: %s... ', i, numElements, name);

    % Sanitize name for struct field
    % 1. Replace non-alphanumeric with underscores
    safeName = regexprep(name, '[^a-zA-Z0-9]', '_');
    % 2. Check if valid start (must be letter)
    if isempty(safeName)
        safeName = 'unnamed_signal';
    elseif ~isletter(safeName(1))
        safeName = ['s' safeName];
    end

    % Use recursive helper to extract data (handles Timeseries, Structs/Buses)
    extracted = process_values(elem.Values);

    if ~isempty(extracted)
        data.(safeName) = extracted;
    else
        fprintf('    Skipped (No valid data found)\n');
    end
end

% --- 4. Save Output ---
[p, n, ~] = fileparts(inputFile);
outputFile = fullfile(p, [n '_py.mat']);

save(outputFile, 'data', '-v7');
fprintf('Successfully converted. Saved to:\n%s\n', outputFile);

end

% --- Helper Functions ---

function out = process_values(val)
% Recursively convert Timeseries or Struct (Bus) to simpler struct
out = [];

if isa(val, 'timeseries')
    % Standard Timeseries
    out.time = val.Time;
    d = val.Data;
    % Flatten 3D arrays (Time, Dim)
    if ndims(d) == 3
        d = permute(d, [3, 1, 2]);
    end
    out.data = d;
    fprintf('OK (timeseries)\n');

elseif isstruct(val)
    % Structure (likely a Bus)
    fprintf('Processing Struct/Bus...\n');
    out = struct();
    fields = fieldnames(val);
    hasData = false;

    for k = 1:numel(fields)
        fn = fields{k};
        subVal = val.(fn);

        fprintf('    -> Field "%s": ', fn);
        res = process_values(subVal);

        if ~isempty(res)
            out.(fn) = res;
            hasData = true;
        else
            fprintf('Skipped\n');
        end
    end

    if ~hasData
        out = [];
    end

elseif isobject(val) && isprop(val, 'Time') && isprop(val, 'Data')
    % Generic Object with Time/Data
    out.time = val.Time;
    out.data = val.Data;
    fprintf('OK (Object with Time/Data)\n');

elseif isstruct(val) && isfield(val, 'time') && isfield(val, 'signals')
    % StructWithTime legacy format
    out.time = val.time;
    out.data = val.signals.values;
    fprintf('OK (StructWithTime)\n');

else
    fprintf('Unrecognized type: %s\n', class(val));
end
end
