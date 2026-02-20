%% Harbour Entry Trajectory Generator
% Generates a trajectory for the "harbourEntry" scenario.
% Shape:
% 1. 5km straight (North)
% 2. Curve Right (Semi-ellipse, 180 deg turn)
%    - Width (Y displacement): 1.65 km
%    - Depth (X displacement): 720 m
% 3. 500m straight (South)
%
% Features:
% - Adjustable ramping of speed at start.
% - Adjustable speed for corner and straight.
% - Tangential heading.

clear; clc;

%% 1. Simulation Parameters
dt = 0.1;              % Sampling time (s)

% Speed Parameters (Adjustable)
v_straight = 2.5;      % Target speed for straight sections (m/s)
v_corner = 2;        % Target speed for corner (m/s)
accel_ramp = 0.05;     % Acceleration for initial ramp (m/s^2) - Slow ramp
accel_change = 0.2;    % Acceleration for changing speed between segments (m/s^2)

% Geometry Parameters
L_straight1 = 5000;    % Length of first straight (m)
L_straight2 = 500;     % Length of second straight (m)

% Curve Parameters
% "1.65km right" -> Total change in Y = 1650m.
% Since it's a 180 turn, this is the diameter in Y.
Ry = 1650 / 2; % 825 m

% "720m forwards and back" -> Max change in X relative to start of curve.
Rx = 720;      % 720 m

% Initial State
x0 = 0;
y0 = 0;
psi0 = 0; % 0 degrees (North)

%% 2. Generate Trajectory Points

% Initialize State
t = 0;
x = x0;
y = y0;
psi = psi0; % radians
v = 0;      % Start still
s = 0;      % Distance traveled

% Buffers for data
% Estimate size: Total distance approx 8000m, avg speed 2m/s -> 2000s. dt=0.1 -> 20000 steps.
% We'll allocate 50000 to be safe.
N_est = 50000;
t_vec = zeros(1, N_est);
x_vec = zeros(1, N_est);
y_vec = zeros(1, N_est);
psi_vec = zeros(1, N_est); % radians
x_dot_vec = zeros(1, N_est);
y_dot_vec = zeros(1, N_est);
psi_dot_vec = zeros(1, N_est); % rad/s

idx = 0; % Counter

% Simulation State
current_segment = 1; % 1: Straight1, 2: Curve, 3: Straight2
finished = false;

% Curve parameter phi (0 to pi)
phi = 0;

fprintf('Generating trajectory...\n');

while ~finished
    % --- 1. Determine Target Speed ---
    v_target = v_straight;

    if current_segment == 1
        % Check distance to corner to decelerate
        dist_to_corner = L_straight1 - x; % Since moving along X
        % Kinematic braking distance: d = (v^2 - v_target^2) / (2*a)
        if v > v_corner
            braking_dist = (v^2 - v_corner^2) / (2 * accel_change);
            if dist_to_corner <= braking_dist + 5 % 5m buffer
                v_target = v_corner;
            end
        end

    elseif current_segment == 2
        v_target = v_corner;

    elseif current_segment == 3
        v_target = v_straight;
    end

    % --- 2. Update Speed (Ramp) ---
    if v < v_target
        v = v + accel_ramp * dt;
        if v > v_target, v = v_target; end
    elseif v > v_target
        v = v - accel_change * dt;
        if v < v_target, v = v_target; end
    end

    % --- 3. Calculate Derivatives and Update Position ---

    if current_segment == 1
        % Straight North
        x_dot = v;
        y_dot = 0;
        psi = 0;
        psi_dot = 0;

        % Update
        x = x + x_dot * dt;
        y = y + y_dot * dt;
        t = t + dt;

        % Transition Check
        if x >= L_straight1
            x = L_straight1; % Clamp
            current_segment = 2;
            phi = 0;
        end

    elseif current_segment == 2
        % Semi-Ellipse
        % Center (Cx, Cy)
        Cx = L_straight1;
        Cy = Ry;
        % Parametric equations (Rotated to start North, turn Right)
        % We want start at (L_straight1, 0) with heading North (0)
        % End at (L_straight1, 2*Ry) with heading South (pi)
        % Max X at L_straight1 + Rx

        % Equations:
        % x = Cx + Rx * sin(phi)
        % y = Cy - Ry * cos(phi)

        % Derivatives w.r.t phi
        dx_dphi = Rx * cos(phi);
        dy_dphi = Ry * sin(phi);

        % Metric ds/dphi
        ds_dphi = sqrt(dx_dphi^2 + dy_dphi^2);

        % dphi/dt = v / (ds/dphi)
        dphi_dt = v / ds_dphi;

        % Derivatives w.r.t time
        x_dot = dx_dphi * dphi_dt;
        y_dot = dy_dphi * dphi_dt;

        % Heading
        psi = atan2(dy_dphi, dx_dphi);

        % Heading Rate (psi_dot)
        % psi = atan2(y', x')
        % dpsi/dphi = (x'y'' - y'x'') / (x'^2 + y'^2)
        % derivatives here are w.r.t phi
        ddx_dphi2 = -Rx * sin(phi);
        ddy_dphi2 = Ry * cos(phi);

        dpsi_dphi = (dx_dphi * ddy_dphi2 - dy_dphi * ddx_dphi2) / (ds_dphi^2);
        psi_dot = dpsi_dphi * dphi_dt;

        % Update
        phi = phi + dphi_dt * dt;
        t = t + dt;

        % Calculate Position from phi
        x = Cx + Rx * sin(phi);
        y = Cy - Ry * cos(phi);

        % Transition Check
        if phi >= pi
            phi = pi;
            current_segment = 3;
            % Clamp position
            x = L_straight1;
            y = 2 * Ry;
            psi = pi;
        end

    elseif current_segment == 3
        % Straight South
        % Moving -X
        x_dot = -v;
        y_dot = 0;
        psi = pi;
        psi_dot = 0;

        % Update
        x = x + x_dot * dt;
        y = y + y_dot * dt;
        t = t + dt;

        % Finish Check
        % End X = L_straight1 - L_straight2
        if x <= (L_straight1 - L_straight2)
            finished = true;
        end
    end

    % Store Data
    idx = idx + 1;
    if idx > length(t_vec)
        % Extend if needed (shouldn't happen with generous preallocation)
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

% Environment (Zero current)
V_c_data = zeros(size(t_vec));
beta_c_data = zeros(size(t_vec));

% Create Timeseries
ts_eta_x       = timeseries(x_vec', t_vec,       'Name', 'eta_x');
ts_eta_y       = timeseries(y_vec', t_vec,       'Name', 'eta_y');
ts_eta_psi     = timeseries(psi_deg', t_vec, 'Name', 'eta_psi'); % Degrees

ts_V_current   = timeseries(V_c_data', t_vec,    'Name', 'V_current');
ts_beta_current= timeseries(beta_c_data', t_vec, 'Name', 'beta_current');

% Derivatives (Feed-forward)
ts_eta_x_dot   = timeseries(x_dot_vec', t_vec,   'Name', 'eta_x_dot');
ts_eta_y_dot   = timeseries(y_dot_vec', t_vec,   'Name', 'eta_y_dot');
ts_eta_psi_dot = timeseries(psi_dot_deg', t_vec, 'Name', 'eta_psi_dot');

%% 4. Create Scenario Object
harbourEntry = Simulink.SimulationData.Dataset;
harbourEntry.Name = 'harbourEntry';

harbourEntry = harbourEntry.addElement(ts_eta_x, 'eta_x');
harbourEntry = harbourEntry.addElement(ts_eta_y, 'eta_y');
harbourEntry = harbourEntry.addElement(ts_eta_psi, 'eta_psi');
harbourEntry = harbourEntry.addElement(ts_V_current, 'V_current');
harbourEntry = harbourEntry.addElement(ts_beta_current, 'beta_current');
% Optional: Add derivatives
harbourEntry = harbourEntry.addElement(ts_eta_x_dot, 'eta_x_dot');
harbourEntry = harbourEntry.addElement(ts_eta_y_dot, 'eta_y_dot');
harbourEntry = harbourEntry.addElement(ts_eta_psi_dot, 'eta_psi_dot');

%% 5. Save
% Save in the same directory as this script
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));
% If running as a script (not function), mfilename might be empty or behave differently depending on how it's called.
% But if called via 'run', it should work.
% Fallback if scriptDir is empty (e.g. if copy-pasted into command window)
if isempty(scriptDir)
    scriptDir = 'Scenarios';
    if ~exist(scriptDir, 'dir')
        mkdir(scriptDir);
    end
end

filename = fullfile(scriptDir, 'harbourEntry.mat');

save(filename, 'harbourEntry', '-v7.3');

fprintf('Scenario saved to %s\n', filename);
fprintf('Total Time: %.2f s\n', t);
