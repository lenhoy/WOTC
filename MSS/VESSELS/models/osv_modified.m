function [xdot,U,M,tau_thr] = osv_modified(vessel, x,ui,Vc,betaVc, tau, skip_alloc)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org)
% [xdot,U,M] = osv(x,ui,Vc,betaVc) returns the speed U in m/s (optionally) 
% and the time derivative xdot of the state vector for an Offshore Supply 
% vessel (OSV). The 6x6 mass matrix M is an optionally output, which can be 
% used for control design.The 6-DOF equations of motion arebased on the 
% nonlinear model of Fossen (2021, Eqs. 6.111-6.116) given by
%   
%   eta_dot = J(eta) * nu
%   nu_dot = nu_c_dot + Minv * ( tau_thr +  tau_drag + tau_crossflow...
%          - (CRB + CA + D) * nu_r - G * eta )
%
% where Minv = inv(MRB + MA) and nu_r = nu - nu_c is the relative 
% velocity. The OSV equipped by two bow tunnel thrusters and two stern 
% azimuth thrusters. 
%
% Inputs: 
%   x: state vector x = [ u v w p q r x y z phi theta psi ]'
%   ui: control inputs ui = [n1, n2, n3, n4, alpha1, alpha2] where n1, n2, 
%       n3 and n4 are the propeller speeds (rps). The last two arguments 
%       alpha1 and alpha2 are the azimuth angles (rad)
%   Vc: OPTIONAL current speed
%   betaVc: OPTIONAL current direction (rad)
%
%   The arguments Vc (m/s) and betaVc (rad) are optional arguments for 
%   ocean currents given by
%
%    v_c = [ Vc * cos(betaVc - psi), Vc * sin( betaVc - psi), 0 ] 
% 
% The generalized thrust vector satisfy (Fossen 2021, Section 11.2.1)
%
%   tau_thr = T_thr(alpha) * K_thr * u_thr
%
% where
%
%   u_thr = [bowThruster 1, bowThruster 2, sternAzimuth 1, sternAzimuth 2]
%   alpha = [alpha1, alpha2]'    
%
% and u_thr(i) = abs(ui(i)) * ui(i) for i = 1...4 is the squared propeller 
% speed. The azimuth angles are defined by alpha = [ui(5), ui(6)]'. The OSV 
% main characteristics are displayed by calling the function OSV without 
% input arguments. The function calls are:
%
%   osv();                                 : Display the OSV main data
%   [~,~,M] = osv()                        : Return the 6x6 mass matrix M
%   [xdot,U] = osv(x,ui,Vc,betaVc,alphaVc) : Return xdot and U, 2-D currents
%   [xdot,U] = osv(x,ui)                   : Return xdot and U, no currents 
% 
% Author:    Thor I. Fossen
% Date:      2024-03-25
% Revisions:
%   2024-04-22: Enhanced compatibility with GNU Octave
%   2024-06-07: Added M as an optional output argument



% Flag for plotting of the surge resitance, linear + quadratic damping
if nargin > 0
    flag = 0;                              
else
    flag = 1;
end

% Initialize state variables
nu = zeros(6,1);
eta = zeros(6,1);

if nargin == 2 || nargin == 0
    Vc = 0;
    betaVc = 0;
end

if nargin == 0                          % Display main ship characteristics
    ui = zeros(6,1);
else
    nu = x(1:6);                        % Generalized velocity vector
    eta = x(7:12);                      % Generalized position vector
end

% Output arguments
M = vessel.M; % Mass matrix
U = sqrt( nu(1)^2 + nu(2)^2 ); % Speed

% Ocean current
v_c = [ Vc * cos(betaVc - eta(6))
        Vc * sin(betaVc - eta(6))
        0                          ];
nu_c = [v_c' zeros(1,3) ]';
nu_c_dot = [-Smtrx(nu(4:6)) * v_c
             zeros(3,1)           ];
nu_r = nu - nu_c;

% Coriolis matrices
[~,CRB] = rbody(vessel.m,vessel.R44,vessel.R55,vessel.R66,nu(4:6),vessel.r_bg'); 
CA  = m2c(vessel.MA, nu);   

% Add linear and quadratic drag in surge using the blending function sigma
[X,Xuu,Xu] = forceSurgeDamping(flag,nu_r(1),vessel.m,vessel.S,vessel.L, ...
    vessel.T1,vessel.rho,vessel.U_max,vessel.thrust_max);

% ---------- ADJUSTMENTS TO REMOVE QUADRATIC DAMPING START ----------------

% set to zero if testing model setup with only linear damping
tau_drag = [ X; zeros(5,1) ];

% Avoid double counting, linear and quadratic damping terms
% comment out if testing model setup with only linear damping
vessel.D(1,1) = 0; % using: X = sigma * Xu * u_r + (1 - sigma) * Xuu * abs(u_r)*u_r

% Add crossflow drag
% set to zero if testing model setup with only linear damping
tau_crossflow = crossFlowDrag(vessel.L,vessel.B,vessel.T,nu_r);

% ---------- ADJUSTMENTS TO REMOVE QUADRATIC DAMPING END ------------------



if skip_alloc
    % Alternatively direct requested thrust straight from controller, skipping alloc. 
    tau_thr = [ tau(1) tau(2) 0 0 0 tau(3) ]';
else
    % Thrust: T = K_max * u, where u = abs(n/n_max) * (n/n_max) is normalized
    u_thr = abs(ui(1:4)) .* ui(1:4);               % Quadratic propeller speed
    alpha = ui(5:6);                               % Azimuth angles
    T_thr = thrConfig( {'T', 'T', alpha(1), alpha(2)}, vessel.l_x, vessel.l_y);
    tau_3dof = T_thr * vessel.K_thr * u_thr;
    tau_thr = [ tau_3dof(1) tau_3dof(2) 0 0 0 tau_3dof(3) ]';
end

% Kinematics
J = eulerang(eta(4),eta(5),eta(6));

% Equations of motion (Fossen 2021, Eqs. 6.111-6.116)
eta_dot = J * nu;
nu_dot = nu_c_dot + ...
    vessel.Minv * ( tau_thr +  tau_drag + tau_crossflow...
    - (CRB + CA + vessel.D) * nu_r - vessel.G * eta);

xdot = [nu_dot; eta_dot];    


%% Print vessel data
if nargin == 0 && nargout == 0

    % Natural frequencies
    w3 = sqrt( vessel.G(3,3) / vessel.M(3,3) );
    w4 = sqrt( vessel.G(4,4) / vessel.M(4,4) );
    w5 = sqrt( vessel.G(5,5) / vessel.M(5,5) );
    T3 = 2 * pi / w3;
    T4 = 2 * pi / w4;    
    T5 = 2 * pi / w5;

    fprintf('\n');
    fprintf('%s\n','-------------------------------------------------------------------------------------');
    fprintf('%s\n','OFFSHORE SUPPLY VESSEL MAIN CHARACTERISTICS');
    fprintf('%s\n','-------------------------------------------------------------------------------------');
    fprintf('%-40s %8.2f m \n', 'Length (L):', vessel.L);
    fprintf('%-40s %8.2f m \n', 'Beam (B):', vessel.B);
    fprintf('%-40s %8.2f m \n', 'Draft (T):', vessel.T);
    fprintf('%-40s %8.2f kg \n', 'Mass (m):', vessel.m);
    fprintf('%-40s %8.2f kg/m^3 \n', 'Density of water (rho):', vessel.rho);
    fprintf('%-40s %8.2f m^3 \n', 'Volume displacement (nabla):', vessel.nabla);
    fprintf('%-40s %8.2f \n', 'Block coefficient (C_b):', vessel.Cb);
    fprintf('%-40s %8.2f \n', 'Waterplane area coefficient (C_w):', vessel.Cw);
    fprintf('%-40s [%2.1f %2.1f %2.1f] m \n', 'Center of gravity (r_bg):',...
        vessel.r_bg(1), vessel.r_bg(2), vessel.r_bg(3));
    fprintf('%-40s %8.2f \n', 'Relative damping ratio in roll:', vessel.zeta4);
    fprintf('%-40s %8.2f \n', 'Relative damping ratio in pitch:', vessel.zeta5);
    fprintf('%-40s %8.2f m \n', 'Transverse metacentric height (GM_T):', vessel.GM_T);
    fprintf('%-40s %8.2f m \n', 'Longitudinal metacentric height (GM_L):', vessel.GM_L);
    fprintf('%-40s %8.2f s \n', 'Natural period in heave (T3):', T3);
    fprintf('%-40s %8.2f s \n', 'Natural period in roll (T4):', T4);
    fprintf('%-40s %8.2f s \n', 'Natural period in pitch (T5):', T5);   
    fprintf('%-40s %8.2f \n', 'Linear surge damping coefficient (Xu):', Xu); 
    fprintf('%-40s %8.2f \n', 'Quadratic drag coefficient (X|u|u):', Xuu); 
   
    vessel.D(1,1) = -Xu;
    matrices = {'Mass matrix: M = MRB + MA', vessel.M;...
        'Linear damping matrix: D', vessel.D; 'Restoring matrix: G', vessel.G};
    
    for k = 1:size(matrices, 1)
        fprintf('%s\n','-------------------------------------------------------------------------------------');
        fprintf('%-40s\n', matrices{k, 1});
        for i = 1:size(matrices{k, 2}, 1)
            for j = 1:size(matrices{k, 2}, 2)
                if matrices{k, 2}(i,j) == 0
                    fprintf('         0 ');
                else
                    fprintf('%10.2e ', matrices{k, 2}(i,j));
                end
            end
            fprintf('\n');
        end
    end

end

end