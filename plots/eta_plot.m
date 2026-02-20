function eta_plot(filePath)
% PLOTSIMRESULTS - Loads data from a simulation .mat file and creates plots.
% If no filePath is provided, it opens a file browser.

% --- 1. Handle Input ---
% If the function is called without a file path, open a dialog to select a file
if nargin == 0
    [fileName, folderPath] = uigetfile('*.mat', 'Select a Simulation Result File');
    if fileName == 0 % User clicked cancel
        disp('No file selected. Plotting cancelled.');
        return;
    end
    filePath = fullfile(folderPath, fileName);
end

% --- 2. Load the Data ---
% The saved file contains a variable named 'out'
load(filePath, 'out');

% --- 3. Extract Data from the Loaded Structure ---
% This assumes the 'out' structure and 'eta' signal exist
try
    % Access the 'eta' element from the logsout Dataset
    eta_element = out.logsout.get('eta');

    % Check if the element was found
    if isempty(eta_element)
        error('Element "eta" not found in logsout.');
    end

    % Extract Time from the Timeseries object
    time = eta_element.Values.Time(:);

    % Access raw data.
    % Error message indicated shape [3 1 90001], which means [Dims x 1 x Time]
    raw_data = eta_element.Values.Data;

    % Squeeze to remove singleton dimensions (e.g. [3 1 N] -> [3 N])
    flat_data = squeeze(raw_data);

    % Check dimensions after squeeze
    if size(flat_data, 2) == length(time)
        % Structure is [Dims x Time]
        % Assuming standard NED: Row 1 = North, Row 2 = East
        eta_x = flat_data(1, :).';
        eta_y = flat_data(2, :).';
    elseif size(flat_data, 1) == length(time)
        % Structure is [Time x Dims]
        eta_x = flat_data(:, 1);
        eta_y = flat_data(:, 2);
    else
        error('Squeezed data dimensions %s do not match Time length %d.', mat2str(size(flat_data)), length(time));
    end

    % Final enforcement of column vectors
    eta_x = eta_x(:);
    eta_y = eta_y(:);

    % Verify consistency
    if length(eta_x) ~= length(time)
        % This should be caught above but double check
        error('Length mismatch after extraction: eta_x=%d, time=%d', length(eta_x), length(time));
    end

    fprintf('Plotting with %d data points.\n', length(time));

catch ME
    error('Could not extract "eta" data: %s. File: %s', ME.message, filePath);
end

% You'll need to load p0 as well. Since p0 wasn't saved in the 'out'
% object, you might need to load it from your parameter file or define it.
% For this example, we'll hardcode it.
p0 = [50; 0]; % Example circle center
p0_x = p0(1);
p0_y = p0(2);

% --- 4. Create the Plot ---
figure;
hold on;
% Scatter expects x, y, size, color.
% Ensure size (25) is scalar, time is vector same size as x,y
scatter(eta_y, eta_x, 25, time, 'filled');
plot(p0_y, p0_x, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
hold off;

% Formatting
title(['Ship Trajectory from: ' filePath], 'Interpreter', 'none');
xlabel('East Position (m)');
ylabel('North Position (m)');
legend('Ship Trajectory', 'Circle Center (p_0)');
grid on;
axis equal;
cb = colorbar;
ylabel(cb, 'Time (s)');

disp(['Plot created successfully from: ' filePath]);

end