% --- Tiled Layout Example ---

% Create some sample data
x = linspace(0, 10, 100);
y1 = sin(x);
y2 = cos(x);
y3 = x.^2;
y4 = randn(size(x));

% 1. Create a figure and a 2x2 tiled layout
figure;
t = tiledlayout(2, 2, 'TileSpacing', 'compact'); % 'compact' reduces space between plots

% Add a main title for the whole dashboard
title(t, 'My Simulation Dashboard');

% --- Plot 1: Top-Left ---
nexttile;
plot(x, y1, 'r-', 'LineWidth', 1.5);
title('Sine Wave');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% --- Plot 2: Top-Right ---
nexttile;
scatter(y1, y2, 10, x, 'filled'); % Use time for color
title('Phase Portrait');
xlabel('sin(x)');
ylabel('cos(x)');
axis equal;

% --- Plot 3: Bottom-Left ---
nexttile;
plot(x, y3, 'g--', 'LineWidth', 1.5);
title('Quadratic Growth');
xlabel('Time (s)');
ylabel('Value');
grid on;

% --- Plot 4: Bottom-Right ---
nexttile;
histogram(y4, 20);
title('Noise Distribution');
xlabel('Value');
ylabel('Frequency');