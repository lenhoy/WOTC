% --- Plotting Script for Ship Position ---

% 1. Extract Data from the Simulation Output
% Assumes 'eta' is a timeseries object inside 'out'
time = out.eta.Time;
eta_x = out.eta.Data(:, 1); % North position
eta_y = out.eta.Data(:, 2); % East position

% Assume p0 (circle center) is defined in the workspace, e.g., p0 = [50; 0];
p0_x = p0(1);
p0_y = p0(2);


% 2. Create the Plot
figure; % Create a new figure window
hold on; % Hold the plot to draw multiple elements

% Use scatter to plot the ship's path, with color indicating time
scatter(eta_y, eta_x, 25, time, 'filled');

% Plot the circle center as a red cross
plot(p0_y, p0_x, 'rx', 'MarkerSize', 12, 'LineWidth', 2);


% 3. Add Labels and Formatting
title('Ship Position Over Time');
xlabel('East Position (m)');
ylabel('North Position (m)');
legend('Ship Trajectory', 'Circle Center (p_0)');
grid on;
axis equal; % Ensures the scaling is the same on both axes
hold off;

% Add a colorbar to show what the colors represent (time in seconds)
cb = colorbar;
ylabel(cb, 'Time (s)');

disp('Plot created successfully.');