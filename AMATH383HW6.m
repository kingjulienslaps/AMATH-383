    %% AMATH 383 HW 6
    
    %Define variables
    sigma = 10;
    b = 8/3;
    r = 28;
    tspan = linspace(0, 1.558652210716175, 1000);
    
    %%Lorenz system.
    lorenz = @(t, X) [sigma * (X(2) - X(1)); 
                       r*X(1) - X(2) - X(1)*X(3); 
                       X(1) * X(2) - b * X(3)];
    %%\sigma (y - x)
    %%rx - y - xz
    %%xy - bz
    
    %Initial conditions
    X0 = [-13.763610682134201, -19.578751942451796, 27];
    
    %Solve using ode45
    [t, X] = ode45(lorenz, tspan, X0);
    
    %Plot
    figure;
    plot3(X(:,1), X(:,2), X(:,3), 'r', 'LineWidth', 1.5);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('Lorenz System Trajectory in Phase Space','Interpreter','latex');
    grid on;
    view(3);
