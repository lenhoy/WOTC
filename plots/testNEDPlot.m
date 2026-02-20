% --- 1. Handle Input ---

[fileName, folderPath] = uigetfile('*.mat', 'Select a Simulation Result File');
if fileName == 0, disp('Plotting cancelled.'); return; end
filePath = fullfile(folderPath, fileName);


% --- 2. Load and Extract Data ---
load(filePath, 'out');
try
    % Extract eta (ship position) trajectory
    eta_signal = out.logsout.get('eta');
    time = eta_signal.Values.Time;
    eta_data = eta_signal.Values.Data; % This is the 3x1xN array
    x = squeeze(eta_data(1, 1, :));
    y = squeeze(eta_data(2, 1, :));
    psi = squeeze(eta_data(3, 1, :));

    % Extract p0 (target position) trajectory
    p0_signal = out.logsout.get('p0');
    p0_data = p0_signal.Values.Data; % This is also a 3x1xN array
    p0_x = p0_data(:, 1);
    p0_y = p0_data(:, 2);
    
catch ME
    error('Could not process signals in %s.\nError: %s', filePath, ME.message);
end


% --- Call the Plotting Function ---
plotNorthEastTrajectory(x, y, psi, time, 1,[0,0],200);