function r = showrate(delta, error2, errorinf, errorgrid)
% This function plots the log-log relationship between delta and the errors
% (L2 norm, L-infinity norm, and grid error) and computes the convergence rate.
% Inputs:
%   delta - Array of delta values
%   error2 - Array of L2 norm errors
%   errorinf - Array of L-infinity norm errors
%   errorgrid - Array of grid errors
% Outputs:
%   r - Array containing the convergence rates for L2 norm, L-infinity norm, and grid error

% Compute the logarithm of delta and errors
ln_x = log(delta);
ln_y1 = log(error2);
ln_y2 = log(errorinf);
ln_y3 = log(errorgrid);

% Plot the log-log relationships
figure;  % Create a new figure
plot(ln_x, ln_y1, 'k-', ln_x, ln_y2, 'k:', ln_x, ln_y3, 'ko');  % Plot with different styles for each error
xlabel('log(\delta)');  % Label for the x-axis
ylabel('log(Error)');  % Label for the y-axis
legend('L2 norm', 'L-infinity norm', 'Grid error');  % Add a legend
title('Log-Log Plot of Error vs. Delta');  % Title for the plot

% Fit a line to the log-log data to compute the convergence rate
a1 = polyfit(ln_x, ln_y1, 1);  % Fit for L2 norm
a2 = polyfit(ln_x, ln_y2, 1);  % Fit for L-infinity norm
a3 = polyfit(ln_x, ln_y3, 1);  % Fit for grid error

% Extract the slope of the fitted lines (which represents the rate)
r = [a1(1), a2(1), a3(1)];

% Display the rates in the command window
disp('Convergence rates:');
disp(['L2 norm rate: ', num2str(r(1))]);
disp(['L-infinity norm rate: ', num2str(r(2))]);
disp(['Grid error rate: ', num2str(r(3))]);
end
