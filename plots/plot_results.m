%% Plot Simulation Results
% Script for generating high-quality plots for the thesis.

%% 1. Load Data
% You can load from SDI or from a saved workspace file
% simData = load_sdi_data(); 

%% 2. Plotting Configuration
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultLineLineWidth', 1.5);

%% 3. Example Plot: Control Inputs
% Assuming we have thruster signals
if exist('simData', 'var') && isfield(simData, 'tau')
    figure(2); clf;
    tau = simData.tau;
    
    subplot(3,1,1);
    plot(tau.Time, tau.Data(:,1));
    ylabel('Surge Force (N)');
    grid on;
    title('Control Inputs');
    
    subplot(3,1,2);
    plot(tau.Time, tau.Data(:,2));
    ylabel('Sway Force (N)');
    grid on;
    
    subplot(3,1,3);
    plot(tau.Time, tau.Data(:,3));
    ylabel('Yaw Moment (Nm)');
    xlabel('Time (s)');
    grid on;
else
    fprintf('Data not loaded or "tau" signal missing.\n');
end
