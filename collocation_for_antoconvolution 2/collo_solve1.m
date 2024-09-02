function Uh = collo_solve1(y_delta, u_exact, N, c, d, alpha, delta, x)
% This function solves a collocation method problem using a combination of
% nonlinear and linear equations.
% Inputs:
%   y_delta - Perturbed data
%   u_exact - Exact solution function
%   N - Number of intervals in the grid
%   c, d - Collocation points parameters (not directly used in this function)
%   alpha - Regularization parameters (array)
%   delta - Perturbation level (noise level)
%   x - A predefined parameter (not directly used in this function)
% Outputs:
%   Uh - Computed approximate solution

h = 1/N;  % Step size
Uh = zeros(N, 1);  % Initialize the solution vector

%%% Compute initial value U0
U0 = chooseu0(u_exact, delta, h, d, x);  % Call a function to choose the initial value

%%%% Solve the nonlinear equation for the first element of Uh
a = c;  
b = 1 - c;  
y1 = y_delta(1);
% Nonlinear equation to solve
f = @(u)(alpha(1)*u + h*a*u*u - alpha(1)*U0 - y1);
u0 = 2;  % 
options = optimoptions('fsolve', 'TolFun', 1e-8, 'MaxIter', 1000);  
usol = fsolve(f, u0, options);  
Uh(1) = usol; 

%%%% Solve the linear equations for the remaining elements of Uh
for n = 1:N-1
    m1 = 0;  % Initialize summation term
    for l = 0:n-1
        % Accumulate the terms based on previous values of Uh
        m1 = m1 + a*Uh(n-l+1)*Uh(l+1) + b*Uh(n-l)*Uh(l+1);
    end
    g = y_delta(n+1) - h*m1 + alpha(n+1)*U0;  % Compute the right-hand side of the linear equation
    Uh(n+1) = g / (alpha(n+1) + h*a*Uh(1));  % Solve for the next Uh value
end
end
