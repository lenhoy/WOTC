%% Pipeline Inspection Trajectory Generator
% Generates a trajectory for the "pipelineInspection" scenario.
% Shape:
% 1. 5km straight (North/X-axis)
% 2. 5km Curve (Radius 2000m)
%    - Moves at constant speed 2.5 m/s throughout.
%
% Coordinates:
% - Starts at (0,0), Heading 0 (North/X-axis).
% - Straight segment moves along +X axis.
% - Curve turns Right (Clockwise).

clear; clc;

%% 1. Simulation Parameters
dt = 0.1;              % Sampling time (s)

% Speed Parameters (Constant)
v_const = 1.5;         % Constant speed (m/s)

% Current Parameters
V_current = 0.0;       % Current speed (m/s)
beta_current_deg = 270;  % Current direction from (degrees)

% Geometry Parameters
L_straight = 2000;     % Length of straight section (m)
L_curve = 5000;        % Length of curve section (m)
R_curve = 2000;        % Radius of curve (m)

% Initial State
x0 = 0;
y0 = 0;
psi0 = 0; % 0 rad (North/X-axis)

%% 2. Generate Trajectory Points

% Initialize State
t = 0;
x = x0;
y = y0;
psi = psi0;
v = v_const;

% Buffers for data
% Total Distance = 10000m. Speed = 2.5. Time = 4000s.
% Steps = 40000. Safety buffer -> 45000.
N_est = 45000;
t_vec = zeros(1, N_est);
x_vec = zeros(1, N_est);
y_vec = zeros(1, N_est);
psi_vec = zeros(1, N_est); % radians
x_dot_vec = zeros(1, N_est);
y_dot_vec = zeros(1, N_est);
psi_dot_vec = zeros(1, N_est); % rad/s

idx = 0;

% Simulation State
current_segment = 1; % 1: Straight, 2: Curve
finished = false;

% Curve parameter phi (0 to L_curve/R_curve)
phi = 0;
max_phi = L_curve / R_curve; % 5000/2000 = 2.5 rad

fprintf('Generating trajectory for Pipeline Inspection...\n');
fprintf('Params: U=%.2f m/s, Vc=%.2f m/s, betac=%.1f deg\n', v_const, V_current, beta_current_deg);

while ~finished
    % Constant Speed
    v = v_const;

    % --- Calculate Derivatives and Update Position ---

    if current_segment == 1
        % Straight North (+X)
        x_dot = v;
        y_dot = 0;
        psi = 0;
        psi_dot = 0;

        % Update
        x = x + x_dot * dt;
        y = y + y_dot * dt;
        t = t + dt;

        % Transition Check
        if x >= L_straight
            x = L_straight; % Clamp to end of straight
            current_segment = 2;
            phi = 0;
        end

    elseif current_segment == 2
        % Turn Right (Clockwise)
        % Center location relative to start of curve (L_straight, 0)
        % Moving +X initially. Turning Right -> Center is at (L_straight, -R)
        Cx = L_straight;
        Cy = -R_curve;

        % Equations for Clockwise turn starting moving +X:
        % x = Cx + R * sin(phi)
        % y = Cy + R * cos(phi)
        % Wait, let's check Start (phi=0):
        % x = L + 0 = L. Correct.
        % y = -R + R = 0. Correct.
        % Check Tangent (dx, dy):
        % dx ~ cos(phi), dy ~ -sin(phi).
        % At phi=0: dx ~ 1, dy ~ 0. Correct (+X direction).

        % Angular velocity
        % ds = R * dphi
        % v = ds/dt = R * dphi/dt
        % dphi/dt = v / R
        dphi_dt = v / R_curve;

        % Derivatives
        dx_dphi = R_curve * cos(phi);
        dy_dphi = -R_curve * sin(phi);

        x_dot = dx_dphi * dphi_dt;
        y_dot = dy_dphi * dphi_dt;

        % Heading
        % psi = atan2(y_dot, x_dot)
        psi = atan2(y_dot, x_dot);

        % Heading Rate
        % psi = -phi (starts at 0, becomes negative)
        % psi_dot = -dphi_dt
        psi_dot = -dphi_dt;

        % Update
        phi = phi + dphi_dt * dt;
        t = t + dt;

        % Position from phi
        x = Cx + R_curve * sin(phi);
        y = Cy + R_curve * cos(phi);

        % Limit Check
        if phi >= max_phi
            finished = true;
        end
    end

    % Store Data
    idx = idx + 1;
    if idx > length(t_vec)
        % Extend
        t_vec = [t_vec, zeros(1, 10000)];
        x_vec = [x_vec, zeros(1, 10000)];
        y_vec = [y_vec, zeros(1, 10000)];
        psi_vec = [psi_vec, zeros(1, 10000)];
        x_dot_vec = [x_dot_vec, zeros(1, 10000)];
        y_dot_vec = [y_dot_vec, zeros(1, 10000)];
        psi_dot_vec = [psi_dot_vec, zeros(1, 10000)];
    end

    t_vec(idx) = t;
    x_vec(idx) = x;
    y_vec(idx) = y;
    psi_vec(idx) = psi;
    x_dot_vec(idx) = x_dot;
    y_dot_vec(idx) = y_dot;
    psi_dot_vec(idx) = psi_dot;

end

% Trim arrays
t_vec = t_vec(1:idx);
x_vec = x_vec(1:idx);
y_vec = y_vec(1:idx);
psi_vec = psi_vec(1:idx);
x_dot_vec = x_dot_vec(1:idx);
y_dot_vec = y_dot_vec(1:idx);
psi_dot_vec = psi_dot_vec(1:idx);


%% 3. Pack Data for Signal Editor
% Convert to degrees for output
psi_deg = rad2deg(psi_vec);
psi_dot_deg = rad2deg(psi_dot_vec);

% Environment
V_c_data = ones(size(t_vec)) * V_current;
beta_c_data = ones(size(t_vec)) * mod((beta_current_deg + 180), 360); % current from -> towards

% Create Timeseries
ts_eta_x       = timeseries(x_vec', t_vec,       'Name', 'eta_x');
ts_eta_y       = timeseries(y_vec', t_vec,       'Name', 'eta_y');
ts_eta_psi     = timeseries(psi_deg', t_vec, 'Name', 'eta_psi'); % Degrees

ts_V_current   = timeseries(V_c_data', t_vec,    'Name', 'V_current');
ts_beta_current= timeseries(beta_c_data', t_vec, 'Name', 'beta_current');

%% 4. Create Scenario Object
% Use a generic name or specific
pipelineInspection = Simulink.SimulationData.Dataset;
pipelineInspection.Name = 'pipelineInspection';

pipelineInspection = pipelineInspection.addElement(ts_eta_x, 'eta_x');
pipelineInspection = pipelineInspection.addElement(ts_eta_y, 'eta_y');
pipelineInspection = pipelineInspection.addElement(ts_eta_psi, 'eta_psi');
pipelineInspection = pipelineInspection.addElement(ts_V_current, 'V_current');
pipelineInspection = pipelineInspection.addElement(ts_beta_current, 'beta_current');

%% 5. Save
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = 'Scenarios';
    if ~exist(scriptDir, 'dir')
        mkdir(scriptDir);
    end
end

% Construct dynamic filename
% Replace decimal points with 'p'
str_U = strrep(sprintf('%.1f', v_const), '.', 'p');
str_Vc = strrep(sprintf('%.1f', V_current), '.', 'p');
str_Beta = sprintf('%d', round(beta_current_deg));

filename_base = sprintf('pipelineInspection_U%s_Vc%s_deg%s.mat', str_U, str_Vc, str_Beta);
filename = fullfile(scriptDir, filename_base);

save(filename, 'pipelineInspection', '-v7.3');

fprintf('Scenario saved to %s\n', filename);
fprintf('Total Time: %.2f s\n', t);
