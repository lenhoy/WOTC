% Plot data from the 10th last simulation run

% 1. Locate the 10th last simulation folder
output_dir = 'output data';
if ~isfolder(output_dir)
    error('Output directory "%s" not found. Run this script from the "Matlab_simulasjon" folder.', output_dir);
end

files = dir(fullfile(output_dir, '20*')); % Filter for timestamped folders starting with 20
dirFlags = [files.isdir];
subDirs = files(dirFlags);

% Sort by name (which is timestamp)
[~, idx] = sort({subDirs.name});
sortedDirs = subDirs(idx);

% Check if we have enough runs
if length(sortedDirs) < 10
    warning('Found only %d runs. Plotting the oldest one available.', length(sortedDirs));
    targetDir = sortedDirs(1);
else
    % Select the 10th last run
    targetDir = sortedDirs(end-9);
end

fprintf('Plotting data from: %s\n', targetDir.name);

% 2. Load the data
matFile = fullfile(output_dir, targetDir.name, [targetDir.name, '.mat']);
if ~isfile(matFile)
    error('Data file not found: %s', matFile);
end

data = load(matFile);

% Handle 'out' structure if present
if isfield(data, 'out')
    simData = data.out;
else
    simData = data;
end

% List available variables
vars = fieldnames(simData);
fprintf('Available variables: %s\n', strjoin(vars, ', '));

% Create figure
figure('Name', ['Simulation: ', targetDir.name], 'Color', 'w');

% 1. Eta (Position)
subplot(2,2,1);
[t_eta, eta] = extract_data(simData, 'eta');
if ~isempty(eta)
    plot(t_eta, eta);
    title('eta (Position)');
    grid on;
    legend('x', 'y', 'psi');
    xlabel('Time (s)');
else
    text(0.5, 0.5, 'eta not found', 'HorizontalAlignment', 'center');
end

% 2. Nu (Velocity)
subplot(2,2,2);
[t_nu, nu] = extract_data(simData, 'nu');
if ~isempty(nu)
    plot(t_nu, nu);
    title('nu (Velocity)');
    grid on;
    legend('u', 'v', 'r');
    xlabel('Time (s)');
else
    text(0.5, 0.5, 'nu not found', 'HorizontalAlignment', 'center');
end

% 3. Tau (Forces)
subplot(2,2,3);
[t_tau, tau] = extract_data(simData, 'tau');
if ~isempty(tau)
    plot(t_tau, tau);
    title('tau (Forces)');
    grid on;
    xlabel('Time (s)');
else
    text(0.5, 0.5, 'tau not found', 'HorizontalAlignment', 'center');
end

% 4. Eta vs Etad (XY Plot)
subplot(2,2,4);
[~, etad] = extract_data(simData, 'etad');

if ~isempty(eta)
    % Plot Actual
    % Assuming eta is [x, y, psi] or similar.
    % Standard: x=North, y=East. Plot East(y) vs North(x).
    if size(eta, 2) >= 2
        plot(eta(:,2), eta(:,1), 'b', 'LineWidth', 1.5); hold on;
        legend_entries = {'Actual'};

        if ~isempty(etad) && size(etad, 2) >= 2
            plot(etad(:,2), etad(:,1), 'r--', 'LineWidth', 1.5);
            legend_entries{end+1} = 'Desired';
        end

        title('XY Plot (North-East)');
        xlabel('East (m)');
        ylabel('North (m)');
        legend(legend_entries);
        grid on;
        axis equal;
    else
        text(0.5, 0.5, 'eta dims < 2', 'HorizontalAlignment', 'center');
    end
else
    text(0.5, 0.5, 'eta not found', 'HorizontalAlignment', 'center');
end


% Local Helper Function
function [t, y] = extract_data(simData, varName)
t = [];
y = [];

if ~isfield(simData, varName)
    return;
end

var = simData.(varName);

if isa(var, 'timeseries')
    t = var.Time;
    y = var.Data;
elseif isstruct(var) && isfield(var, 'time') && isfield(var, 'signals')
    t = var.time;
    y = var.signals.values;
elseif isnumeric(var)
    y = var;
    t = 1:size(y, 1); % Dummy time if not provided
end
end
