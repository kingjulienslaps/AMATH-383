%% AMATH 383 HW4 1E

% Parameters for the exponential growth model
Q0= 1;        % Initial Condition
t = 0:0.01:2;   % Time vector from 0 to 10 with step size 0.1

% Exponential growth model equation
P = 0.75 * Q0 * exp(-2 * t) + Q0 * (-0.75) * exp(-6 * t);
Q = 0.75 * Q0 * exp(-2 * t) + Q0 * (+0.25) * exp(-6 * t);

% Plotting the exponential growth model
figure;
hold on;
h1 = plot(t, P, 'b-', 'LineWidth', 2,'Color','blue');
h2 = plot(t, Q, 'b-', 'LineWidth', 2,'Color','red');

legend([h1, h2], {'$P(t)$, Conc. in bloodstream', '$Q(t)$, Conc. in target muscle'}, ...
       'Location', 'northeast','Interpreter','latex');

xlabel('$t$','Interpreter','latex');
ylabel('Concentration','Interpreter','latex');
title('HW 4 1E Concentration vs Time','Interpreter','latex');
xlim([0, 2]);
ylim([0, 1]);
grid on;
hold off;
%% AMATH 383 HW4 3C
clear
% Define the grid
[P, Q] = meshgrid(-5:0.5:5, -5:0.5:5);  % Adjust step size for clarity

% Define the system (replace with actual dynamics)
% Example: dP/dt = f(P, Q), dQ/dt = g(P, Q)
dP = P -2*Q;   % Example system (adjust as needed)
dQ = -2*P + Q;

% Normalize the vectors
magnitude = sqrt(dP.^2 + dQ.^2);
dP_unit = dP ./ magnitude;
dQ_unit = dQ ./ magnitude;

% Plot the vector field
figure;
quiver(P, Q, dP_unit, dQ_unit, 'b'); % 'b' for blue arrows
hold on;

% Define and plot the trajectory
t = linspace(0, 1, 2); % Adjust time
P_traj = 0.5 * (1) * exp(-t) + 0.5 * (1) * exp(3*t);
Q_traj = 0.5 * (1) * exp(-t) + 0.5 * (-1) * exp(3*t);

% Mark and label the initial condition
plot(1, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % Black dot at (1,0)
text(1.2, 0.2, '(1, 0)', 'FontSize', 12, 'Color', 'k','Interpreter','latex');   % Label with slight offset

plot(P_traj, Q_traj, 'r', 'LineWidth', 2); % 'r' for red trajectory

% Labels and title
xlabel('P','Interpreter','latex');
ylabel('Q','Interpreter','latex');
title('AMATH 383 HW 4 3C. Phase Plane with Vector Field and Trajectory','Interpreter','latex');
axis equal; % To maintain aspect ratio
xlim([-5,5]);
ylim([-5,5]);
grid on;
hold off;

%% AMATH 383 HW 4 3C
% Define the grid
[P, Q] = meshgrid(-5:0.5:5, -5:0.5:5);  % Adjust step size for clarity

% Define the system (replace with actual dynamics)
% Example system: dP/dt = f(P, Q), dQ/dt = g(P, Q)
dP = P-2*Q;   % Example system (adjust as needed)
dQ = -2*P + Q;

% Normalize the vectors
magnitude = sqrt(dP.^2 + dQ.^2);
dP_unit = dP ./ magnitude;
dQ_unit = dQ ./ magnitude;

% Plot the vector field
figure;
h1 = quiver(P, Q, dP_unit, dQ_unit, 'b'); % 'b' for blue arrows
hold on;

% Define and plot the trajectory with initial condition P(0)=1, Q(0)=0
t = linspace(0, 2, 10); % Adjust time range as needed
P_traj = 0.5 * exp(-t) + 0.5 * exp(3*t);
Q_traj = 0.5 * exp(-t) - 0.5 * exp(3*t);

h2 = plot(P_traj, Q_traj, 'r', 'LineWidth', 2); % 'r' for red trajectory

% Mark and label the initial condition
h3 = plot(1, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % Black dot at (1,0)
text(1.2, 0.2, '(1, 0)', 'FontSize', 12, 'Color', 'k');        % Label with slight offset

% Add legend
legend([h1, h2, h3], {'Vector Field', 'Trajectory', 'Initial Condition'}, ...
       'Location', 'northwest','Interpreter','latex');

% Labels and title
xlabel('$P$','Interpreter','latex');
ylabel('$Q$','Interpreter','latex');
title('HW 4 3C Phase Plane with Vector Field and Trajectory','Interpreter','latex');
axis equal; % Maintain aspect ratio

% Set zoom to [-5, 5] in both axes
xlim([-5, 5]);
ylim([-5, 5]);

grid on;
hold off;


%% AMATH 383 HW 4 3E
% Define the grid
[P, Q] = meshgrid(-5:0.5:5, -5:0.5:5);  % Adjust step size for clarity

% Define the system (replace with actual dynamics)
% Example system: dP/dt = f(P, Q), dQ/dt = g(P, Q)
dP = -1*P+2*Q;   % Example sstem (adjust as needed)
dQ = -2*P + Q;

% Normalize the vectors
magnitude = sqrt(dP.^2 + dQ.^2);
dP_unit = dP ./ magnitude;
dQ_unit = dQ ./ magnitude;

% Plot the vector field
figure;
h1 = quiver(P, Q, dP_unit, dQ_unit, 'b'); % 'b' for blue arrows
hold on;

% Define and plot the trajectory with initial condition P(0)=1, Q(0)=0
t = linspace(0,6, 100); % Adjust time range as needed

C1 = 0;
C2 = 2/(sqrt(3));
b = sqrt(3);
a = 0;
%P_traj = C1*(0.5*cos(b*t) - (sqrt(3)*0.5)*sin(b*t)) + C2*(0.5*sin(b*t) + (sqrt(3)*0.5)*cos(b*t));
%Q_traj = C1*cos(b*t)+C2*sin(b*t);
P_traj = C2*(0.5*sin(b*t) + (sqrt(3)*0.5)*cos(b*t));
Q_traj = C2*sin(b*t);

%P_traj = 0.5 * exp(-i*t) + 0.5 * exp(3*t);
%Q_traj = 0.5 * exp(-t) - 0.5 * exp(3*t);

h2 = plot(P_traj, Q_traj, 'r', 'LineWidth', 2); % 'r' for red trajectory

% Mark and label the initial condition
h3 = plot(1, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % Black dot at (1,0)
text(1.2, 0.2, '(1, 0)', 'FontSize', 12, 'Color', 'k','Interpreter','latex');        % Label with slight offset

% Add legend
legend([h1, h2, h3], {'Vector Field', 'Trajectory', 'Initial Condition'}, ...
       'Location', 'northwest','Interpreter','latex');

% Labels and title
xlabel('$P$','Interpreter','latex');
ylabel('$Q$','Interpreter','latex');
title('HW 4 3E Phase Plane with Vector Field and Trajectory','Interpreter','latex');
axis equal; % Maintain aspect ratio

% Set zoom to [-5, 5] in both axes
xlim([-5, 5]);
ylim([-5, 5]);

grid on;
hold off;
