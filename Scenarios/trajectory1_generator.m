%% 1. Simulation Parameters
T_sim = 1000;          % Total simulation time (s)
dt = 0.1;              % Sampling time (s)
t = 0:dt:T_sim;        % Time vector

%% 2. Configure the Ramp
% We want to accelerate from 0 to Target Speed over 't_ramp' seconds.
t_ramp = 180;           % Duration of acceleration phase (s)

% Path parameters
Rx = 200;              % Radius in X (North)
Ry = 150;               % Radius in Y (East)
period = 1000;         % Target time for full loop at cruising speed
w_target = 2*pi / period; % Target angular velocity (rad/s)

%% 3. Calculate Dynamic Angle (Theta)
% Instead of theta = w * t, we calculate theta based on the ramp.
% Phase 1 (t < ramp):  Speed increases linearly (Constant Accel)
% Phase 2 (t >= ramp): Speed is constant (w_target)

theta = zeros(size(t));
theta_dot = zeros(size(t));
theta_ddot = zeros(size(t));

for i = 1:length(t)
    ti = t(i);
    
    if ti < t_ramp
        % --- ACCELERATION PHASE ---
        % Velocity: Linear increase from 0 to w_target
        % v(t) = (w_target / t_ramp) * t
        theta_dot(i) = (w_target / t_ramp) * ti;
        
        % Position: Integral of velocity (1/2 * a * t^2)
        theta(i) = 0.5 * (w_target / t_ramp) * ti^2;
        
        % Acceleration: Constant
        theta_ddot(i) = w_target / t_ramp;
        
    else
        % --- CRUISING PHASE ---
        % Calculate angle at end of ramp to ensure continuity
        theta_ramp_end = 0.5 * w_target * t_ramp;
        
        % Velocity: Constant
        theta_dot(i) = w_target;
        
        % Position: Linear progress from the ramp-end angle
        theta(i) = theta_ramp_end + w_target * (ti - t_ramp);
        
        % Acceleration: Zero
        theta_ddot(i) = 0;
    end
end

%% 4. Calculate Trajectory (Chain Rule)
% Now we map the dynamic Theta into the Oval Geometry
% x = Rx * sin(theta)
% y = Ry * cos(theta) - Ry

% --- POSITION ---
x = Rx * sin(theta);
y = Ry * cos(theta) - Ry; % Shifted to start at [0,0]

% --- VELOCITY (Chain Rule) ---
% dx/dt = dx/dtheta * dtheta/dt
x_dot = (Rx * cos(theta)) .* theta_dot;
y_dot = (-Ry * sin(theta)) .* theta_dot;

% --- ACCELERATION (Product Rule + Chain Rule) ---
% d2x/dt2 = -Rx*sin(theta)*theta_dot^2 + Rx*cos(theta)*theta_ddot
x_ddot = -Rx * sin(theta) .* theta_dot.^2 + Rx * cos(theta) .* theta_ddot;
y_ddot = -Ry * cos(theta) .* theta_dot.^2 - Ry * sin(theta) .* theta_ddot;

%% 5. Heading and Heading Rate
% --- HEADING (psi) ---
% At t=0, velocity is 0, so atan2(0,0) is 0 (North). Perfect.
psi_rad = atan2(y_dot, x_dot);

% Handle the singularity at t=0 if numerical noise occurs (force North)
if x_dot(1) == 0 && y_dot(1) == 0
    psi_rad(1) = 0; 
end

psi_rad = unwrap(psi_rad);
psi_deg = rad2deg(psi_rad);

% --- HEADING RATE (r) ---
% Analytic derivative: (x'y'' - y'x'') / (x'^2 + y'^2)
num = x_dot .* y_ddot - y_dot .* x_ddot;
den = x_dot.^2 + y_dot.^2;

% Handle singularity at start (0/0)
r_rad_s = zeros(size(t));
% Find where velocity is sufficient to calculate curvature
valid_idx = den > 1e-6; 
r_rad_s(valid_idx) = num(valid_idx) ./ den(valid_idx);

r_deg_s = rad2deg(r_rad_s);

%% 6. Pack Data for Signal Editor
% Environment (Zero current)
V_c_data = zeros(size(t));
beta_c_data = zeros(size(t));

% Create Timeseries
ts_eta_x       = timeseries(x', t,       'Name', 'eta_x');
ts_eta_y       = timeseries(y', t,       'Name', 'eta_y');
ts_eta_psi     = timeseries(psi_deg', t, 'Name', 'eta_psi'); % Degrees

ts_V_current   = timeseries(V_c_data', t,    'Name', 'V_current');
ts_beta_current= timeseries(beta_c_data', t, 'Name', 'beta_current');

% Derivatives (Feed-forward)
ts_eta_x_dot   = timeseries(x_dot', t,   'Name', 'eta_x_dot');
ts_eta_y_dot   = timeseries(y_dot', t,   'Name', 'eta_y_dot');
ts_eta_psi_dot = timeseries(r_deg_s', t, 'Name', 'eta_psi_dot');

%% 7. Create Scenario Object
largeOval = Simulink.SimulationData.Dataset;
largeOval.Name = 'largeOval';

largeOval = largeOval.addElement(ts_eta_x, 'eta_x');
largeOval = largeOval.addElement(ts_eta_y, 'eta_y');
largeOval = largeOval.addElement(ts_eta_psi, 'eta_psi');
largeOval = largeOval.addElement(ts_V_current, 'V_current');
largeOval = largeOval.addElement(ts_beta_current, 'beta_current');
% Optional: Add derivatives if you want to inspect them
largeOval = largeOval.addElement(ts_eta_x_dot, 'eta_x_dot');

%% 8. Save
filename = 'Scenarios/largeOval_Scenario.mat';
save(filename, 'largeOval', '-v7.3');

fprintf('Scenario saved to %s\n', filename);
fprintf('Ramp Duration: %d seconds\n', t_ramp);
fprintf('Max Surge Speed: %.2f m/s\n', max(abs(x_dot)));