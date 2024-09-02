function [error2, errorinf] = normcompute(u_exact, Uh, h, c, y)
% This function computes the L2 norm (error2) and L-infinity norm (errorinf)
% between the exact solution u_exact and the computed solution Uh.
% Inputs:
%   u_exact - Exact solution function or piecewise function
%   Uh - Computed solution
%   h - Step size
%   c - Collocation points
%   y - Parameter to select the exact solution or error calculation method
% Outputs:
%   error2 - Computed L2 norm of the error
%   errorinf - Computed L-infinity norm of the error

N = 1/h;  % Number of intervals
error2 = 0;  % Initialize L2 norm error
errorinf = 0;  % Initialize L-infinity norm error
syms s  % Define symbolic variable for integration

for n = 0:N-1
    % Define piecewise exact solution if y == 3
    if y == 3
        if n < 0.5/h
            u_exact = @(s) 0.5;
        elseif (0.5/h <= n) && (n < 0.8/h)
            u_exact = @(s) 0.25;
        else
            u_exact = @(s) 0.75;
        end
    end

    % Case for single collocation point
    if size(c, 2) == 1
        Uhn = Uh(n+1);  % Computed solution at point n
        g = @(s) abs(u_exact(n*h + s*h) - Uhn);  % Absolute difference between exact and computed solutions
        
        % Compute L-infinity norm
        [~, minVal] = fminbnd(@(x) -g(x), 0, 1);  % Find minimum of -g to get maximum of g
        errorinf = max(-minVal, errorinf);  % Update L-infinity error
        
        % Compute L2 norm based on the value of y
        if y == 1
            error2 = error2 + h * ((Uhn - 2)^2 - sin(8*pi*h*n)/(16*h*pi) + sin(8*pi*(h + h*n))/(16*h*pi) ...
                     + (sin(4*pi*h*n)*(Uhn - 2))/(2*h*pi) - (sin(4*pi*(h + h*n))*(Uhn - 2))/(2*h*pi) + 1/2);
        elseif y == 2
            error2 = error2 + h * ((-2*n*h^2 + 2*h)^2/3 - (2*h^2*(-h^2*n^2 + 2*h*n + Uhn - 2))/3 ...
                     + (-h^2*n^2 + 2*h*n + Uhn - 2)^2 + h^4/5 - (h^2*(-2*n*h^2 + 2*h))/2 ...
                     + (-2*n*h^2 + 2*h)*(-h^2*n^2 + 2*h*n + Uhn - 2));
        else
            error2 = error2 + int((Uhn - u_exact(s))^2, s, 0, 1) * h;  % Integrate for L2 norm
        end
    
    % Case for two collocation points
    elseif size(c, 2) == 2
        c1 = c(1); c2 = c(2);  % Collocation points
        Uh1 = Uh(2*n+1); Uh2 = Uh(2*n+2);  % Computed solutions at collocation points
        l1 = @(s) (s - c2) / (c1 - c2);  % Linear basis function 1
        l2 = @(s) (s - c1) / (c2 - c1);  % Linear basis function 2
        f = @(s) abs(Uh1 * l1(s) + Uh2 * l2(s) - u_exact(n*h + s*h));  % Absolute difference between exact and computed solutions
        
        % Compute L-infinity norm
        [~, minVal] = fminbnd(@(x) -f(x), 0, 1);  % Find minimum of -f to get maximum of f
        errorinf = max(-minVal, errorinf);  % Update L-infinity error
        
        % Compute L2 norm based on the value of y
        if y == 1
            error2 = error2 - h * (((sin(4*pi*h*n)*(c1^2 - (Uh1*c2^2)/2 - (Uh2*c1^2)/2 - 2*c1*c2 + c2^2 + (Uh1*c1*c2)/2 + (Uh2*c1*c2)/2)) ...
                    /(h*pi) - (cos(4*pi*h*n)*(Uh1*c1 - Uh1*c2 - Uh2*c1 + Uh2*c2))/(8*h^2*pi^2) + (sin(8*pi*h*n)*(c1^2 - 2*c1*c2 + c2^2)) ...
                    /(16*h*pi))/(c1^2 - 2*c1*c2 + c2^2) - (2*Uh1*c2 - 2*Uh1*c1 + 2*Uh2*c1 - 2*Uh2*c2 - 8*c1*c2 - 4*Uh1*c2^2 - 4*Uh2*c1^2 - Uh1^2*c2 ...
                    - Uh2^2*c1 + 4*c1^2 + 4*c2^2 + (Uh1 - Uh2)^2/3 + Uh1^2*c2^2 + Uh2^2*c1^2 + (c1 - c2)^2/2 + Uh1*Uh2*c1 + Uh1*Uh2*c2 + 4*Uh1*c1*c2 ...
                    + 4*Uh2*c1*c2 + (sin(8*pi*(h + h*n))*(c1^2 - 2*c1*c2 + c2^2))/(16*h*pi) - (cos(4*pi*(h + h*n))*(Uh1*c1 - Uh1*c2 - Uh2*c1 + Uh2*c2)) ...
                    /(8*h^2*pi^2) - (sin(4*pi*(h + h*n))*(Uh1*c1 - Uh1*c2 - Uh2*c1 + Uh2*c2))/(2*h*pi) - 2*Uh1*Uh2*c1*c2 + (sin(4*pi*(h + h*n))*(c1^2 ...
                    - (Uh1*c2^2)/2 - (Uh2*c1^2)/2 - 2*c1*c2 + c2^2 + (Uh1*c1*c2)/2 + (Uh2*c1*c2)/2))/(h*pi))/(c1^2 - 2*c1*c2 + c2^2));
        elseif y == 2
            error2 = error2 + h * ((-2*n*h^2 + 2*h + Uh1/(c1 - c2) - Uh2/(c1 - c2))^2/3 - (h^2*(-2*n*h^2 + 2*h + Uh1/(c1 - c2) - Uh2/(c1 - c2)))/2 ...
                     + (2*h^2*(h^2*n^2 - 2*h*n + (Uh1*c2)/(c1 - c2) - (Uh2*c1)/(c1 - c2) + 2))/3 + h^4/5 - (-2*n*h^2 + 2*h + Uh1/(c1 - c2) - Uh2/(c1 - c2)) ...
                     *(h^2*n^2 - 2*h*n + (Uh1*c2)/(c1 - c2) - (Uh2*c1)/(c1 - c2) + 2) + (h^2*n^2 - 2*h*n + (Uh1*c2)/(c1 - c2) - (Uh2*c1)/(c1 - c2) + 2)^2);
        else
            error2 = error2 + double(int((Uh1*l1(s) + Uh2*l2(s) - u_exact)^2, s, 0, 1) * h);  % Integrate for L2 norm
        end
    end
end

% Final computation of L2 norm
error2 = double(sqrt(error2));
end
