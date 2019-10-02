% Program by Sreevalsan S Menon(sm2hm@mst.edu)


clc                         % Clear command window
clear                       % clear workspace
rng('default')              % For reproducibility
close all                   % Close all figures


%% Hybrid extended Kalman filter

% Dynamics of system
x0=[300000 -20000 0.001]';      % Initial states                             
rho_0=2;                        % Rho value
g=32.2;                         % Gravitational constant
k=20000;                        % air density altitude relation
Q=0;                            % Process noise
R_k=10000;                      % Measurement noise
M=100000;                       % Horizontal Range
a=100000;                       % altitude

[t_m, x_m]=ode45(@(t,x0) system_equations(t,x0), 0:0.5:30,x0);  % Actual system simulation

y=sqrt(M^2+(x_m(:,1)-a).^2)+sqrt(R_k)*randn(size(x_m(:,1)));    % System output
xhat_plus=x0;                                                   % a posteriori estimate for state x
P_plus=[1000000 0 0;0 4000000 0;0 0 10];                        % a posteriori estimate for covariance P
P_minus=P_plus;                                                 % a priori estimate for covariance P
data_EKF=[xhat_plus' xhat_plus' diag(P_plus)'  diag(P_minus)' 0 0 0]; % Saving the data
c=2;                                                            % Variable initialization

for i=0.001:0.001:30                                            % Loop over time
    xhat_dot=[xhat_plus(2); rho_0*exp(-xhat_plus(1)/k)*xhat_plus(2)^2*xhat_plus(3)/2-g;0];% Estimate state update
    xhat_plus=xhat_plus+xhat_dot*0.001;                         % Update a priori estimate state
    % System Dynamics  at nominal value
    A=[0 1 0;-rho_0*exp(-xhat_plus(1)/k)*xhat_plus(2)^2*xhat_plus(3)/(2*k) rho_0*exp(-xhat_plus(1)/k)*xhat_plus(2)*xhat_plus(3) rho_0*exp(-xhat_plus(1)/k)*xhat_plus(2)^2/2;0 0 0];% system matrix
    L=[1 0 0;0 1 0;0 0 1];  % Process covariance
    H_k=[(xhat_plus(1)-a)/sqrt(M^2+(xhat_plus(1)-a)^2) 0 0];    % Output matrix
    M_k=1;                                                      % Measurement ovariance
        
    P_dot=A*P_minus+P_minus*A'+L*Q*L';  % Updating Pdot
    P_minus=P_minus+P_dot*0.001;        % Updating a priori Covariance matrix P
    % Measurement Update    
    if mod(i,0.5)==0                 % Checking for time 
            K_k=P_minus*H_k'*inv(H_k*P_minus*H_k'+M_k*R_k*M_k');                % Updating Kalman gain estimator matriz
            x_hat_plus=xhat_plus+K_k.*(y(c,1)-sqrt(M^2+(xhat_plus(1)-a)^2));    % Updating a posteriori state estimate
            P_plus=(eye(3)-K_k*H_k)*P_minus*(eye(3)-K_k*H_k)'+K_k*M_k*R_k*M_k'*K_k';    % Updating a posteriori covariance matrix
            data_EKF=[data_EKF; x_hat_plus' xhat_plus' diag(P_plus)' diag(P_minus)' K_k'];     % Saving the data
            P_minus=P_plus;                 % Updating a priori Covariance matrix 
            xhat_plus=x_hat_plus;           % Updating a priori state estimate
            c=c+1;                          % Incrementing the variable
    end                             
end

%% Iterated Kalman filter

% Dynamics of system
x0=[300000 -20000 0.001]';  % Initial states                             
rho_0=2;                    % Rho value
g=32.2;                     % Gravitational constant
k=20000;                    % air density altitude relation
Q=0;                        % Process noise
R_k=10000;                  % Measurement noise
M=100000;
a=100000;

xhat_plus=x0;                                                   % a posteriori estimate for state x
P_plus=[1000000 0 0;0 4000000 0;0 0 10];                        % a posteriori estimate for covariance P
P_minus=P_plus;                                                 % a priori estimate for covariance P
data=[xhat_plus' xhat_plus' diag(P_plus)'  diag(P_minus)' 0 0 0]; % Saving the data
c=2;                                                            % Variable initialization

for i=0.001:0.001:30                                            % Loop over time
    xhat_dot=[xhat_plus(2); rho_0*exp(-xhat_plus(1)/k)*xhat_plus(2)^2*xhat_plus(3)/2-g;0];% Estimate state update
    xhat_plus=xhat_plus+xhat_dot*0.001;                         % Update a priori estimate state
    % System Dynamics  at nominal value
    A=[0 1 0;-rho_0*exp(-xhat_plus(1)/k)*xhat_plus(2)^2*xhat_plus(3)/(2*k) rho_0*exp(-xhat_plus(1)/k)*xhat_plus(2)*xhat_plus(3) rho_0*exp(-xhat_plus(1)/k)*xhat_plus(2)^2/2;0 0 0];% system matrix
    L=[1 0 0;0 1 0;0 0 1];  % Process covariance
    H_k=[(xhat_plus(1)-a)/sqrt(M^2+(xhat_plus(1)-a)^2) 0 0];    % Output matrix
    M_k=1;                                                      % Measurement ovariance
        
    P_dot=A*P_minus+P_minus*A'+L*Q*L';  % Updating Pdot
    P_minus=P_minus+P_dot*0.001;        % Updating a priori Covariance matrix P
    % Measurement Update    
    xhatminus=xhat_plus;
    if mod(i,0.5)==0                    % Checking for time 
        for n=1:10   
            H_k=[(xhat_plus(1)-a)/sqrt(M^2+(xhat_plus(1)-a)^2) 0 0];
            K_k=P_minus*H_k'*inv(H_k*P_minus*H_k'+M_k*R_k*M_k');              % Updating Kalman gain estimator matriz
            x_hat_plus=xhat_plus+K_k.*(y(c,1)-sqrt(M^2+(xhat_plus(1)-a)^2));    % Updating a posteriori state estimate
            P_plus=(eye(3)-K_k*H_k)*P_minus*(eye(3)-K_k*H_k)'+K_k*M_k*R_k*M_k'*K_k';    % Updating a posteriori covariance matrix
            xhat_plus=x_hat_plus;     % Updating a priori state estimate
        end
        data=[data; x_hat_plus' xhatminus' diag(P_plus)' diag(P_minus)' K_k'];     % Saving the data
        P_minus=P_plus;                 % Updating a priori Covariance matrix
        c=c+1;                          % Incrementing the variable
    end                             
end


% Plot Altitude error
figure
plot(t_m,abs(x_m(:,1)-data_EKF(:,1)),'--*',t_m,abs(x_m(:,1)-data(:,1)),'-.o', 'LineWidth',2.5,'MarkerSize',5);set(gca,'FontSize',20);set(gca, 'YScale', 'log');
title('Altitude Error','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Altitude','FontSize', 24);
legend('Extended Kalman Filter','Iterated Kalman Filter','FontSize', 24)

% Plot Velocity error
figure
plot(t_m,abs(x_m(:,2)-data_EKF(:,2)),'--*',t_m,abs(x_m(:,2)-data(:,2)),'-.o', 'LineWidth',2.5,'MarkerSize',5);set(gca,'FontSize',20);set(gca, 'YScale', 'log');
title('Velocity Error','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Velocity','FontSize', 24);
legend('Extended Kalman Filter','Iterated Kalman Filter','FontSize', 24)

% Plot Ballistic Coefficient error
figure
plot(t_m,abs(x_m(:,3)-data_EKF(:,3)),'--*',t_m,abs(x_m(:,3)-data(:,3)),'-.o', 'LineWidth',2.5,'MarkerSize',5);set(gca,'FontSize',20);set(gca, 'YScale', 'log');
title('Ballistic Coefficient Error','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Ball. Coeff.','FontSize', 24);
legend('Extended Kalman Filter','Iterated Kalman Filter','FontSize', 24)