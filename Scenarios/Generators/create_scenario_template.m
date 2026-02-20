%% Scenario Generation Template
% This script generates a scenario file (.mat) for the Simulink simulation.
% It defines the path (waypoints or trajectory), environmental forces,
% and other simulation parameters.
%
% Output: A .mat file in the 'Scenarios' folder containing a Simulink.SimulationData.Dataset.

%% 1. Simulation Parameters
T_sim = 1000;          % Total simulation time (s)
dt = 0.1;              % Sampling time (s)
t = 0:dt:T_sim;        % Time vector

%% 2. Path / Trajectory Definition
% Define your trajectory here. You can use mathematical functions or waypoints.

% Example: Circle
R = 100;               % Radius (m)
period = 500;          % Time to complete one circle (s)
w = 2*pi / period;     % Angular velocity (rad/s)

% Position
x = R * sin(w * t);
y = R * cos(w * t);

% Velocity (Derivatives)
x_dot = R * w * cos(w * t);
y_dot = -R * w * sin(w * t);

% Heading (Tangent to path)
psi_rad = atan2(y_dot, x_dot);
psi_deg = rad2deg(unwrap(psi_rad));

% Heading Rate
% For a circle, heading rate is constant w
r_deg_s = rad2deg(w) * ones(size(t));

%% 3. Environmental Forces
% Define current, wind, and waves here.

% Current Speed (m/s)
V_c_data = 0.5 * ones(size(t)); % Constant 0.5 m/s current

% Current Direction (deg)
beta_c_data = 45 * ones(size(t)); % Coming from North-East

%% 4. Pack Data for Simulink
% Create Timeseries objects. Names must match what the Simulink model expects.

ts_eta_x       = timeseries(x', t,       'Name', 'eta_x');
ts_eta_y       = timeseries(y', t,       'Name', 'eta_y');
ts_eta_psi     = timeseries(psi_deg', t, 'Name', 'eta_psi');

ts_V_current   = timeseries(V_c_data', t,    'Name', 'V_current');
ts_beta_current= timeseries(beta_c_data', t, 'Name', 'beta_current');

% Optional: Feed-forward derivatives
ts_eta_x_dot   = timeseries(x_dot', t,   'Name', 'eta_x_dot');
ts_eta_y_dot   = timeseries(y_dot', t,   'Name', 'eta_y_dot');
ts_eta_psi_dot = timeseries(r_deg_s', t, 'Name', 'eta_psi_dot');

%% 5. Create Dataset Object
scenarioData = Simulink.SimulationData.Dataset;
scenarioData.Name = 'CustomScenario'; % Name of the dataset

scenarioData = scenarioData.addElement(ts_eta_x, 'eta_x');
scenarioData = scenarioData.addElement(ts_eta_y, 'eta_y');
scenarioData = scenarioData.addElement(ts_eta_psi, 'eta_psi');
scenarioData = scenarioData.addElement(ts_V_current, 'V_current');
scenarioData = scenarioData.addElement(ts_beta_current, 'beta_current');
scenarioData = scenarioData.addElement(ts_eta_x_dot, 'eta_x_dot');

%% 6. Save to File
% Save in the Scenarios folder (relative to this script or absolute)
% Assuming this script is in Scenarios/Generators/
outputFolder = fileparts(fileparts(mfilename('fullpath'))); % Go up one level to Scenarios
if isempty(outputFolder)
    outputFolder = 'Scenarios'; % Fallback
end

filename = fullfile(outputFolder, 'Custom_Circle_Scenario.mat');
save(filename, 'scenarioData', '-v7.3');

fprintf('Scenario saved to: %s\n', filename);
fprintf('You can now load this scenario in the Simulink model.\n');
