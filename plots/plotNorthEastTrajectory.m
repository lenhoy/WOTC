function plotNorthEastTrajectory(x, y, psi, time, sizeScale, center, viewRange, ax, timeStep)
% plotNorthEastTrajectory Plots a vessel's trajectory with size and view control.
%
% This function visualizes the path of a vessel using patch objects to show
% orientation with varying opacity.
%
% Inputs:
%   x         - A vector of East positions (e.g., in meters).
%   y         - A vector of North positions (e.g., in meters).
%   psi       - A vector of yaw angles in radians.
%   time      - A vector of time values (e.g., in seconds).
%   sizeScale - (Optional) A number from 1 (small) to 10 (large) to
%               control the symbol size. Default is 5.
%   center    - (Optional) A 2-element vector [x, y] for the plot's center.
%   viewRange - (Optional) A number for the view range. The plot will show
%               center +/- viewRange in both axes.
%   ax        - (Optional) Axes handle to plot into.
%   timeStep  - (Optional) Time interval between vessel shapes. Default 1s.

% --- 0. Handle Optional Inputs ---
if nargin < 5 || isempty(sizeScale)
    sizeScale = 5; % Default size setting
end

setCustomView = false;
if nargin >= 7 && ~isempty(center) && ~isempty(viewRange)
    setCustomView = true;
end

% Check if axes handle is provided
if nargin < 8 || isempty(ax)
    figure;
    ax = gca;
end

% Check if timeStep is provided
if nargin < 9 || isempty(timeStep)
    timeStep = 1; % Default: plot every second
end

% --- 1. Basic Trajectory Plot ---
hold(ax, 'on');

plot(ax, x, y, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Vessel Path');

% --- 2. Plot Configuration ---
% Only set title/labels if we created the figure (or if user wants them overwritten)
% For now, we apply them to the axis.
title(ax, 'Vessel Trajectory', 'FontSize', 14);
xlabel(ax, 'East [m]', 'FontSize', 12);
ylabel(ax, 'North [m]', 'FontSize', 12);
grid(ax, 'on');

% --- 3. Define the Ship Symbol Size ---
% If a custom view is set, base the symbol size on that view's range.
% Otherwise, use the full data range.
if setCustomView
    plotDimension = 2 * viewRange;
else
    % Use the data range on the specific axis, or just the data itself
    plotDimension = max(range(x), range(y));
    if plotDimension == 0; plotDimension = 100; end % Fallback for single point
end

% Establish a baseline length relative to the relevant dimension.
% Using 1000 as per your adjustment.
base_length = plotDimension / 1000;

% Multiply the base unit by the user's intuitive sizeScale.
ship_length = base_length * sizeScale;
ship_width = ship_length / 2;

ship_vertices = [
    ship_length/2,  0;              % Bow
    -ship_length/2, ship_width/2;   % Stern, starboard
    -ship_length/2, -ship_width/2;  % Stern, port
    ];

% --- 4. Determine Where to Plot Ship Symbols ---
% We want to avoid plotting too many symbols if time is dense.
plot_times = floor(time(1)):timeStep:floor(time(end));
if isempty(plot_times) && ~isempty(time)
    plot_times = time; % Fallback for very short simulations
end

plot_indices = [];
for t = plot_times
    [~, idx] = min(abs(time - t));
    plot_indices = [plot_indices, idx];
end
plot_indices = unique(plot_indices);

if isempty(plot_indices)
    plot_indices = 1:length(time); % Fallback
end
if plot_indices(end) ~= length(time)
    plot_indices = [plot_indices, length(time)];
end

% --- 5. Draw the Ship Symbols ---
for i = 1:length(plot_indices)
    idx = plot_indices(i);
    current_x = x(idx);
    current_y = y(idx);
    current_psi = psi(idx);

    R = [cos(current_psi), -sin(current_psi); sin(current_psi), cos(current_psi)];
    rotated_vertices = (R * ship_vertices')';
    translated_vertices = rotated_vertices + [current_x, current_y];

    % Fading: 0.05 to 1.0
    if length(plot_indices) > 1
        alpha = 0.05 + 0.95 * ((i - 1) / (length(plot_indices) - 1));
    else
        alpha = 1.0;
    end

    patch(ax, translated_vertices(:,1), translated_vertices(:,2), 'r', ...
        'FaceAlpha', alpha, 'EdgeColor', 'k', 'HandleVisibility', 'off');
end

% --- 6. Set Final View ---
if setCustomView
    xlim(ax, [center(1) - viewRange, center(1) + viewRange]);
    ylim(ax, [center(2) - viewRange, center(2) + viewRange]);
end

axis(ax, 'equal'); % Apply equal axis scaling AFTER setting limits
% hold(ax, 'off'); % Do not turn hold off, let the caller decide
