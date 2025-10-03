%% AMATH 383 HW 5 Question 3h

% Define parameters
beta = 1;
E_values = [1, 1/2, 1/4, 1/8, 1/16];  % Different energy levels

% Define x range
x = linspace(-2, 2, 1000);  % Adjust the range as needed

figure;
hold on;
grid on;

colors = lines(length(E_values));  % Generate distinct colors
legend_entries = [];  % Store legend handles

for i = 1:length(E_values)
    E = E_values(i);
    y_squared = 2*E - x.^2 - (beta/2) * x.^4; % y^2 = 2E - x^2 - 0.5 * beta * x^4
    
    % Ensure valid (real) y values
    valid_idx = y_squared >= 0;
    x_valid = x(valid_idx);
    y_valid = sqrt(y_squared(valid_idx));  % Positive branch
    y_valid_neg = -y_valid;               % Negative branch

    % Plot both branches with a single legend entry per energy level
    h = plot(x_valid, y_valid, 'Color', colors(i, :), 'LineWidth', 1.5); % Positive branch
    plot(x_valid, y_valid_neg, 'Color', colors(i, :), 'LineWidth', 1.5); % Negative branch

    % Store legend entry (only once per energy level)
    legend_entries(i) = h;
end

xlabel('x');
ylabel('y');
title('3h Phase Plane Trajectories for Different Energy Levels','Interpreter','latex');
legend(legend_entries, arrayfun(@(E) sprintf('E = %.4g', E), E_values, 'UniformOutput', false), 'Location', 'best');
hold off;


%% AMATH 383 HW 5 Question 4a.

% Define the domain for P and Q
[P, Q] = meshgrid(linspace(0,2,20), linspace(0,2,20));

% Define the vector field components
dPdt = P - P .* Q;
dQdt = P .* Q - Q;

% Normalize the vectors for readability
magnitude = sqrt(dPdt.^2 + dQdt.^2);
dPdt = dPdt ./ magnitude;
dQdt = dQdt ./ magnitude;

% Plot the vector field
figure;
quiver(P, Q, dPdt, dQdt, 'b', 'LineWidth', 1);
xlabel('P');
ylabel('Q');
title('4a Normalized Vector Field of the System', 'Interpreter','latex');
xlim([0 2]);
ylim([0 2]);
grid on;
axis equal;
%% AMATH 383 HW 5 4b.
% Numerical Solution with ode45

% Define parameters
alpha = 1; beta = 1; gamma = 1; delta = 1;
tspan = linspace(0, 4*pi, 1000); % Time interval with N = 1000 points
y0 = [1; 1/2]; % Initial conditions: P(0) = 1, Q(0) = 1/2

% Define the system of ODEs
dynamics = @(t, y) [y(1) - y(1) * y(2);   % dP/dt = P - PQ
                    y(1) * y(2) - y(2)];  % dQ/dt = QP - Q

% Solve using ode45
[t, Y] = ode45(dynamics, tspan, y0);

% Extract P(t) and Q(t)
P = Y(:,1);
Q = Y(:,2);

% Plot P and Q over time
figure;
plot(t, P, 'b-', 'LineWidth', 1.5); hold on;
plot(t, Q, 'r-', 'LineWidth', 1.5);
xlabel('Time t');
ylabel('P, Q');
title('4b Numerical Solution of the System using ode45','Interpreter','latex');
legend('P(t)', 'Q(t)', 'Location', 'best');
grid on;

%% AMATH 383 HW 5 4c
% 
% - Phase Plane Trajectory with Vector Field

% Define grid for vector field
[P_grid, Q_grid] = meshgrid(linspace(0,2,20), linspace(0,2,20));

% Define the vector field
dPdt_grid = P_grid - P_grid .* Q_grid;
dQdt_grid = P_grid .* Q_grid - Q_grid;

% Normalize vectors for readability
magnitude = sqrt(dPdt_grid.^2 + dQdt_grid.^2);
dPdt_grid = dPdt_grid ./ magnitude;
dQdt_grid = dQdt_grid ./ magnitude;

% Plot the vector field
figure;
quiver(P_grid, Q_grid, dPdt_grid, dQdt_grid, 'b', 'LineWidth', 1);
hold on;

% Plot the phase trajectory from part (b)
plot(P, Q, 'r-', 'LineWidth', 2);

% Labels and formatting
xlabel('P');
ylabel('Q');
title('4c Phase Plane Trajectory with Vector Field','Interpreter','latex');
legend('Vector Field', 'Trajectory', 'Location', 'best');
xlim([0 2]);
ylim([0 2]);
grid on;
axis equal;
hold off

%% Excercise 4e
% Constants
alpha = 1;
beta = 1;
gamma = 1;
delta = 1;
[P_grid, Q_grid] = meshgrid(linspace(0,2,20), linspace(0,2,20));

% Values for E
E_values = [2.01, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3];

% Create a mesh grid for P and Q
%[P, Q] = meshgrid(linspace(0, 2, 50), linspace(0, 2, 50));

% Create a figure for the plot
figure;
hold on;

% Loop through each value of E to plot the trajectories
for E = E_values
    % Equation: delta * P - gamma * log(P) + beta * Q - alpha * log(Q) = E
    % Rearranged to: F(P, Q) = delta * P - gamma * log(P) + beta * Q - alpha * log(Q) - E
    F = delta * P_grid - gamma * log(P_grid) + beta * Q_grid - alpha * log(Q_grid) - E;
    
    % Contour plot for each E value
    contour(P_grid, Q_grid, F, [0, 0], 'LineWidth', 2);
end

% Define the vector field components from the differential equations
dPdt = P_grid .* (alpha - beta * Q_grid);  % dx/dt for P
dQdt = Q_grid .* (delta * P_grid - gamma);  % dy/dt for Q

% Normalize the vectors for readability
magnitude = sqrt(dPdt.^2 + dQdt.^2);
dPdt = dPdt ./ magnitude;
dQdt = dQdt ./ magnitude;

% Plot the vector field
quiver(P_grid, Q_grid, dPdt, dQdt, 'k', 'LineWidth', 1);

% Set plot labels
xlabel('P');
ylabel('Q');
title('4e. Trajectories for Different E Values with Vector Field','Interpreter','latex');
axis tight;
hold off;


%% Exercise 4e
% Constants
alpha = 1;
beta = 1;
gamma = 1;
delta = 1;
[P_grid, Q_grid] = meshgrid(linspace(0,2,20), linspace(0,2,20));

% Values for E
E_values = [2.01, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3];

% Define a colormap for distinct colors
colors = lines(length(E_values)); % Generates a set of distinguishable colors

% Create a figure for the plot
figure;
hold on;

% Loop through each value of E to plot the trajectories
legend_entries = cell(length(E_values), 1); % Preallocate legend entries
contour_handles = gobjects(length(E_values), 1); % Store contour handles

for i = 1:length(E_values)
    E = E_values(i);
    % Equation: delta * P - gamma * log(P) + beta * Q - alpha * log(Q) = E
    % Rearranged to: F(P, Q) = delta * P - gamma * log(P) + beta * Q - alpha * log(Q) - E
    F = delta * P_grid - gamma * log(P_grid) + beta * Q_grid - alpha * log(Q_grid) - E;
    
    % Contour plot for each E value with different colors
    contour_handles(i) = plot(NaN, NaN, 'Color', colors(i, :), 'LineWidth', 2); % Dummy plot for legend
    contour(P_grid, Q_grid, F, [0, 0], 'LineWidth', 2, 'Color', colors(i, :));
    
    % Store legend entry
    legend_entries{i} = sprintf('E = %.2f', E);
end

% Define the vector field components from the differential equations
dPdt = P_grid .* (alpha - beta * Q_grid);  % dx/dt for P
dQdt = Q_grid .* (delta * P_grid - gamma);  % dy/dt for Q

% Normalize the vectors for readability
magnitude = sqrt(dPdt.^2 + dQdt.^2);
dPdt = dPdt ./ magnitude;
dQdt = dQdt ./ magnitude;

% Plot the vector field
quiver(P_grid, Q_grid, dPdt, dQdt, 'k', 'LineWidth', 1);

% Set plot labels
xlabel('P');
ylabel('Q');
title('4e. Trajectories for Different E Values with Vector Field','Interpreter','latex');

% Add legend with correct colors
legend(contour_handles, legend_entries, 'Location', 'best');

axis tight;
hold off;

