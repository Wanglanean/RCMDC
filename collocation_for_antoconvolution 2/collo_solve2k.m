function Uh = collo_solve2k(y_delta, u_exact, N, c, d, alpha, delta, x)
% This function solves a collocation problem with a kernel function using
% both nonlinear and linear equations for a system with two collocation points.
% Inputs:
%   y_delta - Perturbed data
%   u_exact - Exact solution function
%   N - Number of intervals in the grid
%   c - Collocation points parameters (assumed to be a 2-element vector)
%   d - Dimension parameter 
%   alpha - Regularization parameters (array)
%   delta - Perturbation level (noise level)
%   x - A predefined parameter (not directly used in this function)
% Outputs:
%   Uh - Computed approximate solution as a column vector

h = 1/N;  % Step size
m = size(c, 2);  % Number of collocation points (should be 2)
c1 = c(1); c2 = c(2);  % Collocation points
Uh = zeros(N, m);  % Initialize solution matrix

% Compute initial value U0 using the given exact solution
U0 = chooseu0(u_exact, delta, h, d, x);

%%%% Solve the nonlinear equation for the first time step
y1 = y_delta(1); y2 = y_delta(2); 

% Coefficients for the nonlinear equation (computed by AVIE_xishu,m)
a111 = c1 + (-(5*c1^3)/6 + c2*c1^2)/(c1 - c2)^2;
a112 = ((c1^2*(h*c1^2 + 2*c1))/6 - (c1^2*c2*(c1*h + 3))/6)/(c1 - c2)^2;
a121 = ((c1^2*(-h*c1^2 + 2*c1))/6 + (c1^2*c2*(c1*h - 3))/6)/(c1 - c2)^2;
a122 = c1^3/(6*(c1^2 - 2*c1*c2 + c2^2));
a211 = c2^3/(6*(c1 - c2)^2);
a212 = ((c2^2*(-h*c2^2 + 2*c2))/6 + (c1*c2^2*(c2*h - 3))/6)/(c1 - c2)^2;
a221 = ((c2^2*(h*c2^2 + 2*c2))/6 - (c1*c2^2*(c2*h + 3))/6)/(c1 - c2)^2;
a222 = c2 + (-(5*c2^3)/6 + c1*c2^2)/(c1 - c2)^2;

% Define the nonlinear system of equations
f = @(u)nonliearfun(u, h, alpha(1:2), U0, y1, y2, a111, a112, a121, a122, a211, a212, a221, a222);
u0 = [1.5, 1.5]; 
usol = fsolve(f, u0);  
Uh(1, :) = usol;  

%%%% Solve the linear equations for the remaining time steps
for n = 1:N-1
    m1 = 0; m2 = 0;  % Initialize summation terms
    b = zeros(2, 1);  % Initialize the right-hand side vector for linear system
    
    for l = 0:n-1
        % Compute coefficients a111 to a222 and b111 to b222 based on the kernel function
        a111 = (c1*(c1^2 - 6*c1*c2 + 6*c2^2)*(h*n - 2*h*l + 1))/(6*(c1 - c2)^2);
        a112 = ((c1^2*(2*c1 + c1^2*h - 4*c1*h*l + 2*c1*h*n))/6 - (c1^2*c2*(c1*h - 6*h*l + 3*h*n + 3))/6)/(c1 - c2)^2;
        a121 = ((c1^2*(2*c1 - c1^2*h - 4*c1*h*l + 2*c1*h*n))/6 + (c1^2*c2*(c1*h + 6*h*l - 3*h*n - 3))/6)/(c1 - c2)^2;
        a122 = (c1^3*(h*n - 2*h*l + 1))/(6*(c1 - c2)^2);
        a211 = (c2^3*(h*n - 2*h*l + 1))/(6*(c1 - c2)^2);
        a212 = ((c2^2*(2*c2 - c2^2*h - 4*c2*h*l + 2*c2*h*n))/6 + (c1*c2^2*(c2*h + 6*h*l - 3*h*n - 3))/6)/(c1 - c2)^2;
        a221 = ((c2^2*(2*c2 + c2^2*h - 4*c2*h*l + 2*c2*h*n))/6 - (c1*c2^2*(c2*h - 6*h*l + 3*h*n + 3))/6)/(c1 - c2)^2;
        a222 = (c2*(6*c1^2 - 6*c1*c2 + c2^2)*(h*n - 2*h*l + 1))/(6*(c1 - c2)^2);
        
        b111 = ((c1 - 1)*(h + 2*h*l - h*n - 1)*(c1^2 - 6*c1*c2 + 4*c1 + 6*c2^2 - 6*c2 + 1))/(6*(c1 - c2)^2);
        b112 = -(((c1 - 1)^2*(2*c1 - h - 3*c1*h - 2*h*l + h*n + c1^2*h - 4*c1*h*l + 2*c1*h*n + 1))/6 - (c2*(c1 - 1)^2*(c1*h - 4*h - 6*h*l + 3*h*n + 3))/6)/(c1 - c2)^2;
        b121 = (((c1 - 1)^2*(h - 2*c1 + c1*h + 2*h*l - h*n + c1^2*h + 4*c1*h*l - 2*c1*h*n - 1))/6 - (c2*(c1 - 1)^2*(2*h + c1*h + 6*h*l - 3*h*n - 3))/6)/(c1 - c2)^2;
        b122 = ((c1 - 1)^3*(h + 2*h*l - h*n - 1))/(6*(c1 - c2)^2);
        b211 = ((c2 - 1)^3*(h + 2*h*l - h*n - 1))/(6*(c1 - c2)^2);
        b212 = (((c2 - 1)^2*(h - 2*c2 + c2*h + 2*h*l - h*n + c2^2*h + 4*c2*h*l - 2*c2*h*n - 1))/6 - (c1*(c2 - 1)^2*(2*h + c2*h + 6*h*l - 3*h*n - 3))/6)/(c1 - c2)^2;
        b221 = -(((c2 - 1)^2*(2*c2 - h - 3*c2*h - 2*h*l + h*n + c2^2*h - 4*c2*h*l + 2*c2*h*n + 1))/6 - (c1*(c2 - 1)^2*(c2*h - 4*h - 6*h*l + 3*h*n + 3))/6)/(c1 - c2)^2;
        b222 = ((c2 - 1)*(h + 2*h*l - h*n - 1)*(6*c1^2 - 6*c1*c2 - 6*c1 + c2^2 + 4*c2 + 1))/(6*(c1 - c2)^2);
        
        % Accumulate the terms for m1 and m2 based on previous solutions Uh
        m1 = m1 + a111*Uh(n-l+1,1)*Uh(l+1,1) + a112*Uh(n-l+1,1)*Uh(l+1,2) ...
                + a121*Uh(n-l+1,2)*Uh(l+1,1) + a122*Uh(n-l+1,2)*Uh(l+1,2);
        m1 = m1 + b111*Uh(n-l,1)*Uh(l+1,1) + b112*Uh(n-l,1)*Uh(l+1,2) ...
                + b121*Uh(n-l,2)*Uh(l+1,1) + b122*Uh(n-l,2)*Uh(l+1,2);
        
        m2 = m2 + a211*Uh(n-l+1,1)*Uh(l+1,1) + a212*Uh(n-l+1,1)*Uh(l+1,2) ...
                + a221*Uh(n-l+1,2)*Uh(l+1,1) + a222*Uh(n-l+1,2)*Uh(l+1,2);
        m2 = m2 + b211*Uh(n-l,1)*Uh(l+1,1) + b212*Uh(n-l,1)*Uh(l+1,2) ...
                + b221*Uh(n-l,2)*Uh(l+1,1) + b222*Uh(n-l,2)*Uh(l+1,2);
    end
    
    % Construct the linear system A*Uh(n+1,:) = b
    A = zeros(2, 2);  % Initialize matrix A
    
    % Compute coefficients a111 to a222 for the current step
    a111 = -(c1*(h*n - 1)*(c1^2 - 6*c1*c2 + 6*c2^2))/(6*(c1 - c2)^2);
    a112 = ((c1^2*(2*c1 + c1^2*h - 2*c1*h*n))/6 - (c1^2*c2*(c1*h - 3*h*n + 3))/6)/(c1 - c2)^2;
    a121 = -((c1^2*(c1^2*h - 2*c1 + 2*c1*h*n))/6 - (c1^2*c2*(c1*h + 3*h*n - 3))/6)/(c1 - c2)^2;
    a122 = -(c1^3*(h*n - 1))/(6*(c1 - c2)^2);
    a211 = -(c2^3*(h*n - 1))/(6*(c1 - c2)^2);
    a212 = -((c2^2*(c2^2*h - 2*c2 + 2*c2*h*n))/6 - (c1*c2^2*(c2*h + 3*h*n - 3))/6)/(c1 - c2)^2;
    a221 = ((c2^2*(2*c2 + c2^2*h - 2*c2*h*n))/6 - (c1*c2^2*(c2*h - 3*h*n + 3))/6)/(c1 - c2)^2;
    a222 = -(c2*(h*n - 1)*(6*c1^2 - 6*c1*c2 + c2^2))/(6*(c1 - c2)^2);
    
    % Fill the matrix A
    A(1, 1) = alpha(m*n+1) + h*(a111*Uh(1,1) + a121*Uh(1,2));
    A(1, 2) = h*(a112*Uh(1,1) + a122*Uh(1,2));
    A(2, 1) = h*(a211*Uh(1,1) + a221*Uh(1,2));
    A(2, 2) = alpha(m*n+2) + h*(a212*Uh(1,1) + a222*Uh(1,2));
    
    % Compute the right-hand side vector b
    b(1) = y_delta(m*n+1) - h*m1 + alpha(n*m+1)*U0;
    b(2) = y_delta(m*n+2) - h*m2 + alpha(n*m+2)*U0;
    
    % Solve the linear system for Uh at the next time step
    Uh(n+1, :) = (A\b)';
end

% Reshape Uh to a column vector for output
Uh = reshape(Uh', m*N, 1);
end
