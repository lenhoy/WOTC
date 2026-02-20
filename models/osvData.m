%% OSV Ship model parameters
%global vessel;

vessel.L = 83;                 % Length (m)
vessel.B = 18;                 % Beam (m)
vessel.T = 5;                  % Draft (m)
vessel.rho = 1025;             % Density of water (kg/m3)
vessel.Cb = 0.65;              % Block coefficient: Cb = nabla / (L * B * T)
vessel.S = vessel.L * vessel.B + 2 * vessel.T * vessel.B; % Wetted suwrface, box approximation

% Thrust: T = K_max * abs(n/n_max) * (n/n_max) = K_thr * abs(n) * n
vessel.K_max = [300e3 300e3 655e3 655e3]';% Max propeller thrust (N)
vessel.n_max = [140 140 150 150]'; % Max propeller speed (rpm)
vessel.K_thr = diag(vessel.K_max./vessel.n_max.^2);% Thruster coefficient matrix
vessel.l_x = [37, 35, -vessel.L/2, -vessel.L/2]; % Thruster x-coordinates
vessel.l_y = [0, 0, -7, 7]; % Thruster y-coordinates

vessel.thrust_max = vessel.K_max(3)+vessel.K_max(4); % max thrust in the surge direction (N)
vessel.U_max = 7.7; % Max cruise speed (m/s) corresponding to max thrust

vessel.nabla = vessel.Cb * vessel.L * vessel.B * vessel.T; % Volume displacement(m3)
vessel.m = vessel.rho * vessel.nabla; % Mass (kg)
vessel.r_bg = [-4.5 0 -1.2]'; % Location of the CG with respect to the CO

vessel.Cw = 0.8; % Waterplane area coefficient: Cw = Awp/(L * B)
vessel.Awp = vessel.Cw * vessel.B * vessel.L; % Waterplane area
vessel.KB = (1/3) * (5*vessel.T/2 - vessel.nabla/vessel.Awp); % Eq. (4.38)
vessel.k_munro_smith = (6 * vessel.Cw^3) / ...
    ( (1+vessel.Cw) * (1+2*vessel.Cw)); % Eq. (4.37)
vessel.r_bb = [-4.5 0 vessel.T-vessel.KB]'; % Location of the CB with respect to the CO
vessel.BG = vessel.r_bb(3) - vessel.r_bg(3); % Vertical distance between CG and CB

vessel.I_T = vessel.k_munro_smith * ...
    (vessel.B^3 * vessel.L) / 12; % Transverse moment of inertia
vessel.I_L = 0.7 * (vessel.L^3 * vessel.B) / 12;% Longitudinal moment of inertia
vessel.BM_T = vessel.I_T / vessel.nabla;
vessel.BM_L = vessel.I_L / vessel.nabla;
vessel.GM_T = vessel.BM_T - vessel.BG; % Should be larger than 0.5 m
vessel.GM_L = vessel.BM_L - vessel.BG;

% G matrix
vessel.LCF = -0.5;% x-distance from the CO to the center of Awp
vessel.r_bp = [0 0 0]'; % Compute G in the CO
vessel.G = Gmtrx(vessel.nabla,vessel.Awp,vessel.GM_T,vessel.GM_L,vessel.LCF,vessel.r_bp);

% Rigid-body mass matrix MRB
vessel.R44 = 0.35 * vessel.B; % Radius of gyration in roll, see Eq.(4.77)-(4.78)
vessel.R55 = 0.25 * vessel.L; % Radius of gyration in pitch
vessel.R66 = 0.25 * vessel.L; % Radius of gyration in yaw
vessel.MRB = rbody(vessel.m,vessel.R44,vessel.R55,vessel.R66,[0,0,0]',vessel.r_bg); 

% The added mass matrix MA is derived from the frequency-dependent potential
% coefficients using a look-alike supply vessel in the MSS toolbox. The data
% is stored in the structure <vessel>.
%   load supply             % Check data by typing vessel
%   disp(vessel.main)
%   vessel = computeManeuveringModel(vessel, 1, 7, [3, 1.2, 3.3], 1);
vessel.MA = 1e9 * [0.0006   0    0    0    0    0
     0    0.0020         0    0.0031         0   -0.0091
     0         0    0.0083         0    0.0907         0
     0    0.0031         0    0.0748         0   -0.1127
     0         0    0.0907         0    3.9875         0
     0   -0.0091         0   -0.1127         0    1.2416];

% Mass matrix including hydrodynamic added mass
vessel.M = vessel.MRB + vessel.MA;
vessel.Minv = invQR(vessel.M);

% D matrix
vessel.T1 = 100;               % Time constants for linear damping (s)
vessel.T2 = 100;
vessel.T6 = 1;
vessel.zeta4 = 0.15;           % Relative damping ratio in roll
vessel.zeta5 = 0.3;            % Relative damping ratio in pitch
vessel.D = Dmtrx([vessel.T1, vessel.T2, vessel.T6], ...
    [vessel.zeta4,vessel.zeta5],vessel.MRB,vessel.MA,vessel.G);

%% Setting up a bus to get data into simulink
busInfo = Simulink.Bus.createObject(vessel);
