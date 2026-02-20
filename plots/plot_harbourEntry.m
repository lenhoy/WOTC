%% Plot Harbour Entry Trajectory
% Visualizes the harbourEntry scenario and performs automated verification.

clear; clc; close all;

%% 1. Load Data
% Determine path to scenario file relative to this script
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = 'plots'; % Fallback
end
% Assuming folder structure: Root/plots/this_script.m and Root/Scenarios/file.mat
% So we need to go up one level from plots, then into Scenarios
projectRoot = fileparts(scriptDir);
scenarioFile = fullfile(projectRoot, 'Scenarios', 'harbourEntry.mat');

if ~exist(scenarioFile, 'file')
    error('Scenario file not found: %s', scenarioFile);
end

fprintf('Loading %s...\n', scenarioFile);
load(scenarioFile); % Loads 'harbourEntry' dataset

% Extract Timeseries
ts_x = harbourEntry.getElement('eta_x');
ts_y = harbourEntry.getElement('eta_y');
ts_psi = harbourEntry.getElement('eta_psi');
ts_x_dot = harbourEntry.getElement('eta_x_dot');
ts_y_dot = harbourEntry.getElement('eta_y_dot');

t = ts_x.Time;
x = ts_x.Data;
y = ts_y.Data;
psi_deg = ts_psi.Data;
x_dot = ts_x_dot.Data;
y_dot = ts_y_dot.Data;

% Calculate Speed
U = sqrt(x_dot.^2 + y_dot.^2);

%% 2. Automated Verification
fprintf('\n--- Automated Verification ---\n');
passed = true;

% Metrics
max_x = max(x);
max_y = max(y);
end_x = x(end);
end_y = y(end);
end_psi = psi_deg(end);

% Expected Values (Approximate)
% 5km straight -> Turn starts at x=5000.
% Turn depth 720m -> Max X = 5000 + 720 = 5720.
% Turn width 1.65km -> Max Y = 1650.
% Return 500m -> End X = 5000 - 500 = 4500.
% End Heading -> 180 degrees.

tol_pos = 10; % meters
tol_ang = 1;  % degrees

% Check Max X
if abs(max_x - 5720) < tol_pos
    fprintf('[PASS] Max X: %.2f m (Expected ~5720)\n', max_x);
else
    fprintf('[FAIL] Max X: %.2f m (Expected ~5720)\n', max_x);
    passed = false;
end

% Check Max Y
if abs(max_y - 1650) < tol_pos
    fprintf('[PASS] Max Y: %.2f m (Expected ~1650)\n', max_y);
else
    fprintf('[FAIL] Max Y: %.2f m (Expected ~1650)\n', max_y);
    passed = false;
end

% Check End X
if abs(end_x - 4500) < tol_pos
    fprintf('[PASS] End X: %.2f m (Expected ~4500)\n', end_x);
else
    fprintf('[FAIL] End X: %.2f m (Expected ~4500)\n', end_x);
    passed = false;
end

% Check End Heading
% Normalize to [0, 360) or similar for comparison if needed, but 180 is expected.
if abs(abs(end_psi) - 180) < tol_ang
    fprintf('[PASS] End Heading: %.2f deg (Expected 180)\n', end_psi);
else
    fprintf('[FAIL] End Heading: %.2f deg (Expected 180)\n', end_psi);
    passed = false;
end

if passed
    fprintf('>>> VERIFICATION SUCCESSFUL <<<\n');
else
    fprintf('>>> VERIFICATION FAILED <<<\n');
end

%% 3. Visualization

fig = figure('Name', 'Harbour Entry Trajectory', 'Color', 'w', 'Position', [100, 100, 1200, 800]);

% --- Main Plot: Trajectory ---
ax1 = subplot(2, 2, [1, 3]); % Left side
hold(ax1, 'on');
grid(ax1, 'on');
axis(ax1, 'equal');
xlabel(ax1, 'East (m)'); % Y is East usually, but here X is North, Y is East?
% Script said: x=North, y=East.
xlabel(ax1, 'Y (East) [m]');
ylabel(ax1, 'X (North) [m]');
title(ax1, 'Trajectory (Color = Time)');

% Plot colored line (Surface hack for gradient line)
% z = zeros(size(x));
% surface([y, y], [x, x], [z, z], [t, t], ...
%         'FaceColor', 'no', 'EdgeColor', 'interp', 'LineWidth', 2, 'Parent', ax1);
% Use scatter for simplicity or patch if needed, but surface is standard for gradient line.
scatter(ax1, y, x, 10, t, 'filled'); % Simple scatter with color
colormap(ax1, 'parula');
c = colorbar(ax1);
c.Label.String = 'Time (s)';

% Quiver for heading (Decimated)
decim = 500; % Plot every Nth point
quiver(ax1, y(1:decim:end), x(1:decim:end), ...
    sin(deg2rad(psi_deg(1:decim:end))), cos(deg2rad(psi_deg(1:decim:end))), ...
    0.5, 'k', 'LineWidth', 1);

% Current Position Marker
hMarker = plot(ax1, y(1), x(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% --- Subplot 1: Speed ---
ax2 = subplot(2, 2, 2);
plot(ax2, t, U, 'b-', 'LineWidth', 1.5);
grid(ax2, 'on');
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Speed (m/s)');
title(ax2, 'Speed Profile');
hLineU = xline(ax2, t(1), 'r--', 'LineWidth', 1.5);

% --- Subplot 2: Heading ---
ax3 = subplot(2, 2, 4);
plot(ax3, t, psi_deg, 'k-', 'LineWidth', 1.5);
grid(ax3, 'on');
xlabel(ax3, 'Time (s)');
ylabel(ax3, 'Heading (deg)');
title(ax3, 'Heading Profile');
hLinePsi = xline(ax3, t(1), 'r--', 'LineWidth', 1.5);

%% 4. Interactivity
% Slider
hSlider = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Units', 'normalized', 'Position', [0.2, 0.02, 0.6, 0.03], ...
    'Min', 1, 'Max', length(t), 'Value', 1, ...
    'Callback', @(src, event) updatePlot(src, t, x, y, hMarker, hLineU, hLinePsi));

% Add text label for slider
uicontrol('Parent', fig, 'Style', 'text', 'String', 'Time Slider', ...
    'Units', 'normalized', 'Position', [0.45, 0.05, 0.1, 0.03], ...
    'BackgroundColor', 'w');

%% Callback Function
function updatePlot(slider, t, x, y, hMarker, hLineU, hLinePsi)
idx = round(slider.Value);

% Update Marker
hMarker.XData = y(idx);
hMarker.YData = x(idx);

% Update Vertical Lines
hLineU.Value = t(idx);
hLinePsi.Value = t(idx);

drawnow;
end
