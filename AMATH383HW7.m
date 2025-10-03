%% AMATH 383 HW 7

clc; clear; close all;

% Define parameters
A = 3/2;
B = -1;
N = 10;

% Define function
f = @(p) p.* (A + B* p.^2);

% Generate cobweb diagrams for p0 = 1 and p0 = -1
p0 = [1, -1];
cobweb_plot(f, 1, N, A, B);
cobweb_plot(f, -1, N, A, B);

% Cobweb plot function
function cobweb_plot(f, p0, n_iter, A, B)
    figure; hold on;
    fplot(f, [-1.5, 1.5], 'b', 'LineWidth', 1.5); % Plot f(p)
    fplot(@(p) p, [-1.5, 1.5], 'k--', 'LineWidth', 1.2); % Plot p_n+1 = p_n
    xlabel('p_n'); ylabel('p_{n+1}');
    title(['Cobweb Diagram for p_0 = ', num2str(p0)]);
    grid on;

    % Generate cobweb steps
    p_n = p0;
    for i = 1:n_iter
        p_next = f(p_n);
        plot([p_n, p_n], [p_n, p_next], 'r', 'LineWidth', 1.2); % Vertical line
        plot([p_n, p_next], [p_next, p_next], 'r', 'LineWidth', 1.2); % Horizontal line
        p_n = p_next;
    end
    %plot(-sqrt((1-A)/B),-sqrt((1-A)/B),'o')
    %plot(sqrt((1-A)/B),sqrt((1-A)/B),'o')

    hold off;
end



%% Excercise 5

