% Program by Sreevalsan S Menon(sm2hm@mst.edu)

clc                         % Clear command window
clear                       % clear workspace
rng('default')              % For reproducibility
close all                   % Close all figures


%% Continous-Time Extended Kalman filter 
% Dynamics of system
x0=[100000 -6000 2000]';    % Initial states                             
rho_0=0.0034;               % Rho value
g=32.2;                     % Gravitational constant
k=22000;                    % air density altitude relation
Q=0;                        % Process noise
R=100;                      % Measurement noise
[t, x]=ode45(@(t,x0) system_equations(t,x0), 0:0.004:16,x0); % simulating the actual system

% Simulating system for 100 runs
for sim=1:100                       % Simulating 100 program runs
    y=x+sqrt(R)*randn(size(x));     % Finding output 
    x_hat=[100010 -6100 2500]';     % Initial state estimate
    P=[500 0 0;0 20000 0;0 0 250000]; % Initial Covariance matrix
    data=[ x_hat' diag(P)' 0 0 0];  % Saving data
    
    % Finding the estimate
    for i=1:size(t,1)               % Loop over time
        % System Dynamics  at nominal value
        A=[0 1 0;-rho_0*exp(-x_hat(1)/k)*x_hat(2)^2/(2*k*x_hat(3)) rho_0*exp(-x_hat(1)/k)*x_hat(2)/x_hat(3) -rho_0*exp(-x_hat(1)/k)*x_hat(2)^2/(2*x_hat(3)^2);0 0 0];% System matrix
        L=[1 0 0;0 1 0;0 0 1];      % Process covariance
        C=[1 0 0];                  % Output matrix
        M=1;                        % Measurement covariance 

        Q_tilda=L*Q*L';             % Process noise
        R_tilda=M*R*M';             % Measurement noise

        P_dot=A*P+P*A'+Q_tilda-P*C'*inv(R_tilda)*C*P; % Updation of covariance matrix Pdot
        P=P+P_dot*0.004;            % Updation of covariance matrix P
        K=P*C'*inv(R_tilda);        % Updating Kalman estimator gain
        xhat_dot=[x_hat(2); rho_0*exp(-x_hat(1)/k)*x_hat(2)^2/(2*x_hat(3)) - g;0]+K*(y(i,1)-x_hat(1));% Finding state estimate
        x_hat=x_hat+xhat_dot*0.004; % Finding state estimate
        data=[data; x_hat' diag(P)' K']; % Saving the data

    end
    datasim(:,:,sim)=data;          % Saving the data
    clear data                      % Clear variable
end
data_extended=mean(datasim,3);      % Finding mean value over runs
clear datasim                       % Clear variable

%% Hybrid extended Kalman filter

% Dynamics of system
x0=[100000 -6000 2000]';    % Initial states                             
rho_0=0.0034;               % Rho value
g=32.2;                     % Gravitational constant
k=22000;                    % air density altitude relation
Q=0;                        % Process noise
R_k=100;                    % Measurement noise
[t_m, x_m]=ode45(@(t,x0) system_equations(t,x0), 0:0.5:16,x0);% Actual system simulation

% Simulating system for 100 runs
for sim=1:100                           % Simulating 100 program runs
    
    y=x_m+sqrt(R_k)*randn(size(x_m));   % System output
    xhat_plus=[100010 -6100 2500]';     % a posteriori estimate for state x
    P_plus=[500 0 0;0 20000 0;0 0 250000];% % a posteriori estimate for covariance P
    P_minus=P_plus;                     % a priori estimate for covariance P
    xhat_minus=xhat_plus;               % a priori estimate for state x
    data=[xhat_plus' xhat_minus' diag(P_plus)'  diag(P_minus)' 0 0 0]; % Saving the data
    c=2;                                % Variable initialization
% Finding the estimate
    for i=2:size(t,1)                   % Loop over time
        xhat_dot=[xhat_minus(2); rho_0*exp(-xhat_minus(1)/k)*xhat_minus(2)^2/(2*xhat_minus(3)) - g;0];% Estimate state update
        xhat_minus=xhat_minus+xhat_dot*0.004; % Update a priori estimate state
        % System Dynamics  at nominal value
        A=[0 1 0;-rho_0*exp(-xhat_minus(1)/k)*xhat_minus(2)^2/(2*k*xhat_minus(3)) rho_0*exp(-xhat_minus(1)/k)*xhat_minus(2)/xhat_minus(3) -rho_0*exp(-xhat_minus(1)/k)*xhat_minus(2)^2/(2*xhat_minus(3)^2);0 0 0];% system matrix
        L=[1 0 0;0 1 0;0 0 1];  % Process covariance
        H_k=[1 0 0];            % Output matrix
        M_k=1;                  % Measurement ovariance
        
        P_dot=A*P_minus+P_minus*A'+L*Q*L';  % Updating Pdot
        P_minus=P_minus+P_dot*0.004;        % Updating a priori Covariance matrix P
        
        if mod(t(i),0.5)==0                 % Checking for time 
            K_k=P_minus*H_k'*inv(H_k*P_minus*H_k'+M_k*R_k*M_k');    % Updating Kalman gain estimator matriz
            x_hat_plus=xhat_minus+K_k.*(y(c,1)-H_k*xhat_minus);     % Updating a posteriori state estimate
            P_plus=(eye(3)-K_k*H_k)*P_minus*(eye(3)-K_k*H_k)'+K_k*M_k*R_k*M_k'*K_k';    % Updating a posteriori covariance matrix
            data=[data; x_hat_plus' xhat_minus' diag(P_plus)' diag(P_minus)' K_k'];     % Saving the data
            P_minus=P_plus;                                                             % Updating a priori Covariance matrix 
            xhat_minus=x_hat_plus;          % Updating a priori state estimate
            c=c+1;                          % Incrementing the variable
        end                             
    end
    datasim(:,:,sim)=data;                  % Saving the data
end 
data_hybrid=mean(datasim,3);                % Finding mean for 100 simulations

% Plot Altitude
figure
plot(t(1:100:end),abs(x(1:100:end,1)-data_extended(1:100:end,1)),'--*',t_m,abs(x_m(:,1)-data_hybrid(:,1)),'-.o', 'LineWidth',1.5,'MarkerSize',5);set(gca,'FontSize',20);
title('Altitude','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Altitude','FontSize', 24);
legend('Continous EKF','Hybrid EKF','FontSize', 24)

% Plot Velocity 
figure
plot(t(1:100:end),abs(x(1:100:end,2)-data_extended(1:100:end,2)),'--*',t_m,abs(x_m(:,2)-data_hybrid(:,2)),'-.o', 'LineWidth',1.5,'MarkerSize',5);set(gca,'FontSize',20);
title('Velocity','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Velocity','FontSize', 24);
legend('Continous EKF','Hybrid EKF','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(t(3001:100:end),abs(x(3001:100:end,2)-data_extended(3001:100:end,2)),'--*',t_m(25:end),abs(x_m(25:end,2)-data_hybrid(25:end,2)),'-.o', 'LineWidth',1.5,'MarkerSize',5);set(gca,'FontSize',10)

% Plot 1/Ballistic Coefficient
figure
plot(t(1:100:end),abs(x(1:100:end,3)-data_extended(1:100:end,3)),'--*',t_m,abs(x_m(:,3)-data_hybrid(:,3)),'-.o', 'LineWidth',1.5,'MarkerSize',5);set(gca,'FontSize',20);
title('1/Ballistic Coefficient','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('1/Ball. Coeff.','FontSize', 24);
legend('Continous EKF','Hybrid EKF','FontSize', 24)

%% Find simulation estimation error

err_simulation=[ std(x(1:100:end,1)-data_extended(1:100:end,1)) std(x_m(:,1)-data_hybrid(:,1));std(x(1:100:end,2)-data_extended(1:100:end,2)) std(x_m(:,2)-data_hybrid(:,2));std(x(1:100:end,3)-data_extended(1:100:end,3)) std(x_m(:,3)-data_hybrid(:,3))]; 
fprintf('Standard deviation of estimated error for altitude from simulation for Continous is %d : and Hybrid is : %d feet\n' ,err_simulation(1,1),err_simulation(1,2))
fprintf('Standard deviation of estimated error for Velocity from simulation for Continous is %d : and Hybrid is : %d feet/s\n' ,err_simulation(2,1),err_simulation(2,2))
fprintf('Standard deviation of estimated error for 1/Ballistic coefficient from simulation for Continous is %d : and Hybrid is : %d \n' ,err_simulation(3,1),err_simulation(3,2))

%% Plot a priori and a posteriori  for covariance matrix
ti(3:2:66)=1:32;ti(2:2:65)=1:32;
P_11(1:2:65)=data_hybrid(:,7);P_11(2:2:65)=data_hybrid(2:end,10);
figure
plot(ti,P_11,0:32,data_hybrid(:,7),'*',1:32,data_hybrid(2:end,10),'o','LineWidth',1.5,'MarkerSize',5)
title('Covariance Matrix - a priori and a posteriori altitude estimation error variance(Hybrid)','FontSize', 24)
xlabel('Time steps','FontSize', 24);ylabel('P_{11}','FontSize', 24);
legend('P_{11} ','P_{11}^+','P_{11}^-','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(ti(floor(end/2):end),P_11(floor(end/2):end),16:32,data_hybrid(end/2:end,7),'*',16:32,data_hybrid(end/2:end,10),'o','LineWidth',1.5,'MarkerSize',5)

ti(3:2:66)=1:32;ti(2:2:65)=1:32;
P_22(1:2:65)=data_hybrid(:,8);P_22(2:2:65)=data_hybrid(2:end,11);
figure
plot(ti,P_22,0:32,data_hybrid(:,8),'*',1:32,data_hybrid(2:end,11),'o','LineWidth',1.5,'MarkerSize',5)
title('Covariance Matrix - a priori and a posteriori velocity estimation error variance(Hybrid)','FontSize', 24)
xlabel('Time steps','FontSize', 24);ylabel('P_{22}','FontSize', 24);
legend('P_{22} ','P_{22}^+','P_{22}^-','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(ti(floor(end/2):end),P_22(floor(end/2):end),16:32,data_hybrid(end/2:end,8),'*',16:32,data_hybrid(end/2:end,11),'o','LineWidth',1.5,'MarkerSize',5)


ti(3:2:66)=1:32;ti(2:2:65)=1:32;
P_22(1:2:65)=data_hybrid(:,9);P_22(2:2:65)=data_hybrid(2:end,12);
figure
plot(ti,P_22,0:32,data_hybrid(:,9),'*',1:32,data_hybrid(2:end,12),'o','LineWidth',1.5,'MarkerSize',5)
title('Covariance Matrix - a priori and a posteriori 1/Ballistic coefficient estimation error variance(Hybrid)','FontSize', 24)
xlabel('Time steps','FontSize', 24);ylabel('P_{33}','FontSize', 24);
legend('P_{33} ','P_{33}^+','P_{33}^-','FontSize', 24)
