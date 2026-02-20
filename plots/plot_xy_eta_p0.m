function plot_xy_eta_p0(filePath)
% PLOT_XY_ETA_P0 - Loads and plots simulation data for two trajectories.

% --- 1. Handle Input ---
if nargin == 0
    [fileName, folderPath] = uigetfile('*.mat', 'Select a Simulation Result File');
    if fileName == 0, disp('Plotting cancelled.'); return; end
    filePath = fullfile(folderPath, fileName);
end

% --- 1.5 Configuration ---
plotP0 = false; % Toggle to enable/disable plotting of target p0

% Vessel Plot Settings
shipSizeScale = 10;       % Size of the vessel triangle (default 1-10)
shipPlotInterval = 60;    % Time interval between vessel shapes in seconds (e.g. 60 = once per minute)

% --- 2. Load and Extract Data ---
load(filePath, 'out');
try
    % Extract eta (ship position) trajectory
    eta_signal = out.logsout.get('eta');
    time = eta_signal.Values.Time;
    eta_data = eta_signal.Values.Data; % This is the 3x1xN array
    eta_x = squeeze(eta_data(1, 1, :));
    eta_y = squeeze(eta_data(2, 1, :));
    eta_psi = squeeze(eta_data(3, 1, :)); % Extract Yaw (Psi)

    % Extract p0 (target position) trajectory safely
    p0_x = []; p0_y = [];
    if plotP0
        p0_signal = out.logsout.getElement('p0'); % Use getElement for safety if available, or try/catch
        if isempty(p0_signal)
            p0_signal = out.logsout.get('p0'); % Fallback to get
        end

        if ~isempty(p0_signal)
            p0_data = p0_signal.Values.Data;
            % Handle dimensions: could be Nx2/Nx3 or 3x1xN (if timeseries)
            if ndims(p0_data) == 3
                p0_x = squeeze(p0_data(1, 1, :));
                p0_y = squeeze(p0_data(2, 1, :));
            else
                p0_x = p0_data(:, 1);
                p0_y = p0_data(:, 2);
            end
        else
            disp('Warning: Signal "p0" not found in logsout.');
        end
    end

catch ME
    error('Could not process signals in %s.\nError: %s', filePath, ME.message);
end

% --- 3. Downsample the Data ---
% DEFINE HOW MANY SAMPLES TO SKIP. A larger number means fewer points.
plotInterval = 10;

% Select every Nth point from the full datasets
% p0 data is decimated for the scatter plot
if ~isempty(p0_x)
    p0_x_decimated = p0_x(1:plotInterval:end);
    p0_y_decimated = p0_y(1:plotInterval:end);
end


% --- 4. Create the Plot ---

figure;
ax = gca;
hold(ax, 'on');

% Plot the ship's trajectory (eta) using the specialized function to current axis
% We use the decimated data for the symbols to avoid overcrowding, or full data for path?
% plotNorthEastTrajectory handles its own time-based sampling for symbols,
% so passing full decimated data is fine.
plotNorthEastTrajectory(eta_x, eta_y, eta_psi, time, shipSizeScale, [], [], ax, shipPlotInterval);

% Plot the target's trajectory (p0) with filled squares
if plotP0 && ~isempty(p0_x)
    scatter(ax, p0_x_decimated, p0_y_decimated, 30, 'k', 'x', 'DisplayName', 'Circle Center (p_0)');
end


% --- 5. Formatting ---
[~, fileName, ~] = fileparts(filePath);
title(ax, ['Ship Position in NED with Circle Center: ' fileName], 'Interpreter', 'none');
xlabel(ax, 'East Position (m)');
ylabel(ax, 'North Position (m)');
legend(ax, 'show');
grid(ax, 'on');
axis(ax, 'equal');
% colormap('jet');
% cb = colorbar;
% ylabel(cb, 'Time (s)');

disp(['Plot created successfully from: ' filePath]);

hold off;

end