% Program by Sreevalsan S Menon(sm2hm@mst.edu)

% Method - Recursive least square

clc                                                     % Clear command window
clear                                                   % Clear workspace
close all                                               % Close figures

year=1946:1956;                                         % Year
time=[1:11]';                                           % Year time step
production=[66.6;84.9;88.6; 78.0;96.8; 105.2;93.2; 111.6;88.3;117.0;115.2]; % Production

% Find Linear fit
xhat = [0;0];                                           % Initialize xhat
Pk = 1e6*eye(2);                                        % Initialize P0
Rk=10;                                                  % Initialize Rk
for k=1:size(time,1)                                    % Recursive loop
    Hk=[1 k];                                           % H Matrix
    Kk=Pk*Hk'*inv(Hk*Pk*Hk'+Rk);                        % finding Kk matrix
    xhat=xhat+Kk*(production(k)-Hk*xhat);               % Update xhat
    Pk=(eye(2)-Kk*Hk)*Pk*(eye(2)-Kk*Hk)'+Kk*Rk*Kk';     % Update Pk
end                                                     % Loop terminates
linear=xhat;                                            % Saving linear coefficients

% Find quadratic fit
xhat = [0;0;0];                                         % Initialize xhat
Pk = 1e6*eye(3);                                        % Initialize P0
Rk=10;                                                  % Initialize Rk
for k=1:size(time,1)                                    % Recursive loop
    Hk=[1 k k^2];                                       % H Matrix
    Kk=Pk*Hk'*inv(Hk*Pk*Hk'+Rk);                        % finding Kk matrix
    xhat=xhat+Kk*(production(k)-Hk*xhat);               % Update xhat
    Pk=(eye(3)-Kk*Hk)*Pk*(eye(3)-Kk*Hk)'+Kk*Rk*Kk';     % Update Pk
end                                                     % Loop terminates
quadratic=xhat;                                         % Saving quadratic coefficients

% Find Cubic fit
xhat = [0;0;0;0];                                       % Initialize xhat
Pk = 1e6*eye(4);                                        % Initialize P0
Rk=10;                                                  % Initialize Rk
for k=1:size(time,1)                                    % Recursive loop
    Hk=[1 k k^2 k^3];                                   % H Matrix
    Kk=Pk*Hk'*inv(Hk*Pk*Hk'+Rk);                        % finding Kk matrix
    xhat=xhat+Kk*(production(k)-Hk*xhat);               % Update xhat
    Pk=(eye(4)-Kk*Hk)*Pk*(eye(4)-Kk*Hk)'+Kk*Rk*Kk';     % Update Pk
end                                                     % Loop terminates
cubic=xhat;                                             % Saving Cubic coefficients

% Find quartic fit
xhat = [0;0;0;0;0];                                     % Initialize xhat
Pk = 1e6*eye(5);                                        % Initialize P0
Rk=10;                                                  % Initialize Rk
for k=1:size(time,1)                                    % Recursive loop
    Hk=[1 k k^2 k^3 k^4];                               % H Matrix
    Kk=Pk*Hk'*inv(Hk*Pk*Hk'+Rk);                        % finding Kk matrix
    xhat=xhat+Kk*(production(k)-Hk*xhat);               % Update xhat
    Pk=(eye(5)-Kk*Hk)*Pk*(eye(5)-Kk*Hk)'+Kk*Rk*Kk';     % Update Pk
end                                                     % Loop terminates
quartic=xhat;                                           % Saving quartic coefficients

% Predicting using fit
linear_pred = linear'*[ones(1,11); time'];
quadratic_pred = quadratic'*[ones(1,11); time'; time'.^2];
cubic_pred= cubic'*[ones(1,11); time'; time'.^2; time'.^3];
quartic_pred=quartic'*[ones(1,11); time'; time'.^2; time'.^3; time'.^4];

% Plotting fits
figure
plot(year,production,'*',year,linear_pred)
xlabel('Year');ylabel('Production in million Tons')
title('Linear Curve Fit - Recursive least square')
figure
plot(year,production,'*',year,quadratic_pred)
xlabel('Year');ylabel('Production in million Tons')
title('Quadratic Curve Fit - Recursive least square')
figure
plot(year,production,'*',year,cubic_pred)
xlabel('Year');ylabel('Production in million Tons')
title('Cubic Curve Fit - Recursive least square')
figure
plot(year,production,'*',year,quartic_pred)
xlabel('Year');ylabel('Production in million Tons')
title('Quartic Curve Fit - Recursive least square')

% Finding RMS Error
RMSE_linear = sqrt(mean((production' - linear_pred).^2));
RMSE_quadratic = sqrt(mean((production' - quadratic_pred).^2));
RMSE_cubic = sqrt(mean((production' - cubic_pred).^2));
RMSE_quartic = sqrt(mean((production' - quartic_pred).^2));

% Print out RMS Error
fprintf('The RMS error for linear curve fit using recursive least square is: %d \n',RMSE_linear)
fprintf('The RMS error for quadratic curve fit using recursive least square is: %d \n',RMSE_quadratic)
fprintf('The RMS error for cubic curve fit using recursive least square is: %d \n',RMSE_cubic)
fprintf('The RMS error for quartic curve fit using recursive least square is: %d \n\n',RMSE_quartic)

% Prediction for 1957
linear_pred_year = linear'*[1; 12];
quadratic_pred_year = quadratic'*[1; 12 ;12^2];
cubic_pred_year = cubic'*[1; 12; 12^2; 12^3];
quartic_pred_year =quartic'*[1 ;12; 12^2 ;12^3 ;12^4];

% Print out prediction for 1957
fprintf('Prediction for 1957 using linear curve fit using recursive least square is: %d \n',linear_pred_year)
fprintf('Prediction for 1957 using quadratic curve fit using recursive least square is: %d \n',quadratic_pred_year)
fprintf('Prediction for 1957 using cubic curve fit using recursive least square is: %d \n',cubic_pred_year)
fprintf('Prediction for 1957 using quartic curve fit using recursive least square is: %d \n',quartic_pred_year)