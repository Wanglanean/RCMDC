clc; clear;
h = 0.001;  % Step size for grid discretization

%%%% Selection of collocation points parameters (Gaussian collocation nodes)
c = 1/2;  % Alternative collocation node
%c = [(3-sqrt(3))/6, (3+sqrt(3))/6];  % Gaussian collocation nodes
m = size(c, 2);  % Number of collocation nodes
d = m;  % Dimension of collocation nodes

%%%% Mesh discretization
N = round(1/h);  % Number of intervals in the grid
h = 1/N;  % Recalculate step size based on number of intervals
T = zeros(N*m, 1);  % Initialize time grid
for n = 0:N-1
    for j = 1:m
        T(m*n + j) = n*h + c(j)*h;  % Compute time grid points based on collocation nodes
    end
end

%%%% Regularization parameters
%%% Variable regularization parameters
alpha_min = 0.3; alpha_max = 0.4;  % Define the range for regularization parameters
alpha_increase = linspace(alpha_min, alpha_max, N*m);  % Generate a linearly spaced array of alpha values
load('randorderm1.mat', 'alpha_randorder');
%load('randorderm2.mat', 'alpha_randorder');
alpha_adp = alpha_increase(alpha_randorder);  % Randomly permute the alpha values for adaptation

%%% Fixed regularization parameters
nopt = 20;  % Number of fixed regularization parameters
alpha_opt = zeros(N*m, nopt);  % Initialize matrix to hold fixed alpha values
for i = 1:nopt
    alpha_fix = alpha_increase(round(linspace(1, N*m, nopt)));  % Select alpha values from the generated array
    alpha_opt(:, i) = linspace(alpha_fix(i), alpha_fix(i), N*m);  % Replicate the fixed alpha value across the grid
end

%%%% Define exact solution
%%% Continuous function 
u_exact = @(t) 2 + cos(4*pi*t);  % Cosine-based exact solution
%u_exact = @(t) t.^2 - 2*t + 2;  % Quadratic exact solution
%%% Piecewise function
%u_exact = @(t) (t >= 0 & t <= 0.5) * 0.5 + (t > 0.5 & t <= 0.8) * 0.25 + (t > 0.8 & t <= 1) * 0.75;  % Piecewise exact solution
u = u_exact(T);  

%%%% Kernel function 
k = @(t,s) 1;  % Constant kernel function
%k = @(t,s) t - 2*s + 1;  % Non-constant kernel function

%%%% Load perturbed data y_delta
y_exaxt = computey(u_exact, k, T);
load('y_delta_Ex1m1delta1.mat', 'y_delta');
delta = max(abs(y_delta - y_exaxt));  % Compute the maximum perturbation (noise level)
%%%% Collocation method solution
x = 1; y = 1;  % Some predefined parameters for collocation

%%% Solution with variable regularization parameter
% Solve using variable regularization
Uh_adp = collo_solve1(y_delta, u_exact, N, c, d, alpha_adp, delta, x); 
% Plot exact solution and approximate solution
figure(1)
plot(T, u, "k-", T, Uh_adp, "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12)  
%print('Figure/Ex3m2delta2', '-dpng', '-r600');
% Compute error norms for variable regularization
[erroradp2, erroradpinf] = normcompute(u_exact, Uh_adp, h, c, y);  
erroradpgrid = norm(u - Uh_adp, inf);  

%%% Solution with fixed regularization parameter
Uh_opt = zeros(N*m, nopt); 
error_opt_2 = zeros(1, nopt); error_opt_inf = zeros(1, nopt); error_opt_grid = zeros(1, nopt); 
for i = 1:nopt
    % Solve using fixed regularization
    Uh_opt(:, i) = collo_solve1(y_delta, u_exact, N, c, d, alpha_opt(:, i), delta, x); 
    % Compute error norms for fixed regularization
    [error_opt_2(i), error_opt_inf(i)] = normcompute(u_exact, Uh_opt(:, i), h, c, y);  
    error_opt_grid(i) = norm(u - Uh_opt(:, i), inf);
end
[min_erroropt2, min_erroropt2_index] = min(error_opt_2);  % Find minimum L2 error and its index
[min_erroroptinf, min_erroroptinf_index] = min(error_opt_inf);  
[min_erroroptgrid, min_erroroptgrid_index] = min(error_opt_grid);
% Plot exact solution and best fixed regularization solution
figure(2)
plot(T, u, "k-", T, Uh_opt(:, min_erroropt2_index), "k--") 
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12) 
%print('Figure/Ex3m2delta2alpha_fix', '-dpng', '-r600'); 
