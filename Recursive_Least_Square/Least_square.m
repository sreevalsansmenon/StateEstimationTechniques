% Program by Sreevalsan S Menon(sm2hm@mst.edu)

% Method least square

clc                                                     % Clear command window
clear                                                   % Clear workspace
close all                                               % Close figures

year=1946:1956;                                         % Year
time=[1:11]';                                           % Year time step
production=[66.6;84.9;88.6; 78.0;96.8; 105.2;93.2; 111.6;88.3;117.0;115.2]; % Production

H=[ones(11,1) time];                                    % Linear H Matrix
linear = inv(H'*H)*H'*production                       % Finding linear coefficients
H=[ones(11,1) time time.^2];                            % Quadratic H Matrix
quadratic = inv(H'*H)*H'*production                    % Finding quadratic coefficients
H=[ones(11,1) time time.^2 time.^3];                    % Cubic H Matrix
cubic = inv(H'*H)*H'*production                        % Finding cubic coefficients
H=[ones(11,1) time time.^2 time.^3 time.^4];            % Quartic H Matrix
quartic = inv(H'*H)*H'*production                      % Finding quartic coefficients

% Predicting using fit
linear_pred = linear'*[ones(1,11); time'];
quadratic_pred = quadratic'*[ones(1,11); time'; time'.^2];
cubic_pred= cubic'*[ones(1,11); time'; time'.^2; time'.^3];
quartic_pred=quartic'*[ones(1,11); time'; time'.^2; time'.^3; time'.^4];

% Plotting fits
figure
plot(year,production,'*',year,linear_pred)
xlabel('Year');ylabel('Production in million Tons')
title('Linear Curve Fit')
figure
plot(year,production,'*',year,quadratic_pred)
xlabel('Year');ylabel('Production in million Tons')
title('Quadratic Curve Fit- Least Square')
figure
plot(year,production,'*',year,cubic_pred)
xlabel('Year');ylabel('Production in million Tons')
title('Cubic Curve Fit - Least Squar')
figure
plot(year,production,'*',year,quartic_pred)
xlabel('Year');ylabel('Production in million Tons')
title('Quartic Curve Fit - Least Square')

% Finding RMS Error
RMSE_linear = sqrt(mean((production' - linear_pred).^2));
RMSE_quadratic = sqrt(mean((production' - quadratic_pred).^2));
RMSE_cubic = sqrt(mean((production' - cubic_pred).^2));
RMSE_quartic = sqrt(mean((production' - quartic_pred).^2));

% Print out RMS Error
fprintf('The RMS error for linear curve fit using least square is: %d \n',RMSE_linear)
fprintf('The RMS error for quadratic curve fit using least square is: %d \n',RMSE_quadratic)
fprintf('The RMS error for cubic curve fit using least square  is: %d \n',RMSE_cubic)
fprintf('The RMS error for quartic curve fit using least square  is: %d \n\n',RMSE_quartic)

% Prediction for 1957
linear_pred_year = linear'*[1; 12];
quadratic_pred_year = quadratic'*[1; 12 ;12^2];
cubic_pred_year = cubic'*[1; 12; 12^2; 12^3];
quartic_pred_year =quartic'*[1 ;12; 12^2 ;12^3 ;12^4];

% Print out prediction for 1957
fprintf('Prediction for 1957 using linear curve fit using least square is: %d \n',linear_pred_year)
fprintf('Prediction for 1957 using quadratic curve fit using least square is: %d \n',quadratic_pred_year)
fprintf('Prediction for 1957 using cubic curve fit using least square is: %d \n',cubic_pred_year)
fprintf('Prediction for 1957 using quartic curve fit using least square is: %d \n',quartic_pred_year)