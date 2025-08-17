function y = computey(u, k, T)
% This function computes the values of y based on the integral of a given integrand.
% Inputs:
%   u - Function handle representing the exact solution or a function to be integrated
%   k - Function handle representing the kernel function
%   T - Array of time points at which to evaluate the integrals
% Outputs:
%   y - Array of computed values based on the integrals at each time point in T

y = zeros(length(T), 1);  % Initialize the output array y

for i = 1:length(T)
    t = T(i);  % Current time point
    % Define the integrand as a function of s, where the integrand depends on t, u, and k
    integrand = @(s) k(t, s) .* u(t - s) .* u(s);  
    % Compute the integral of the integrand from 0 to t
    y(i) = integral(integrand, 0, t);  
end
end
