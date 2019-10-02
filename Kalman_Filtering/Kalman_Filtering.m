% Program by Sreevalsan S Menon(sm2hm@mst.edu)

clc                         % Clear command window
clear                       % clear workspace
rng('default')              % For reproducibility
close all                   % Close all figures

%% Dynamics of system

F_k_minus_1=[0.5 2;0 1];    % System matrix                             
G_k_minus_1=[0 0;0 0];      % Input matrix
u_k_minus_1=[0;0];          % Input
H_k=[1 0];                  % Output matrix
Q_k=[0 0;0 10];             % Covariance matrix process
R_k=10;                     % Covariance matrix measurement 

%% Kalman filter initialization

xhat_k_plus=[600;200];      % a posteriori estimate for state x
p_k_plus=[500 0;0 200];     % a posteriori estimate for covariance P
x_k=[650;250];              % initial state
data=[x_k' xhat_k_plus' 0 0 diag(p_k_plus)' 0 0 0 0]; % Vector to save data over loop

%% Kalmna filtering

for k=1:10                                                          % Loop to execute k timepoints
    
    x_k=F_k_minus_1*x_k+[0 ; sqrt(10)*randn];                       % Actual system dynamics 
    y_k(k)=H_k*x_k+sqrt(10)*randn;                                  % Actual system output
    p_k_minus=F_k_minus_1*p_k_plus*F_k_minus_1'+Q_k;                % Updating a priori estimate for covariance P
    K_k=p_k_minus*H_k'*inv(H_k*p_k_minus*H_k'+R_k);                 % Kalman estimator gain matrix
    xhat_k_minus=F_k_minus_1*xhat_k_plus+G_k_minus_1*u_k_minus_1;   % Updating a priori estimate for state x
    y_hat_k(k)=H_k*xhat_k_minus;                                    % Output estimate
    xhat_k_plus=xhat_k_minus+K_k*(y_k(k)-H_k*xhat_k_minus);         % Updating a posteriori estimate for state x
    p_k_plus=(eye(2)-K_k*H_k)*p_k_minus*(eye(2)-K_k*H_k)'+K_k*R_k*K_k';% Updating a posteriori estimate for covariance P
    data=[data; x_k' xhat_k_plus' xhat_k_minus' diag(p_k_plus)' diag(p_k_minus)' K_k'];%  Saving the data in each iteration
    
end

%% Plot True vs Estimated Population
figure
plot(1:k,y_hat_k,'--*',1:k,y_k,'-.o', 'LineWidth',1.5,'MarkerSize',5);set(gca,'FontSize',20);
title('True vs Estimated Population','FontSize', 24)
xlabel('Time Steps','FontSize', 24);ylabel('Population','FontSize', 24);
legend('Estimated Population','True Population','FontSize', 24)

%% Plot True vs Estimated Food Supply
figure
plot(1:k,data(1:end-1,2),'--*',1:k,data(1:end-1,4),'-.o', 'LineWidth',1.5,'MarkerSize',5);set(gca,'FontSize',20);
title('True vs Estimated Food Supply','FontSize', 24)
xlabel('Time Steps','FontSize', 24);ylabel('Food Supply','FontSize', 24);
legend('Estimated Food Supply','True Food Supply','FontSize', 24)

%% Plot Standard Deviation of Estimation Error
figure
plot(1:k,sqrt(data(1:end-1,7)),'--*',1:k,sqrt(data(1:end-1,8)),'-.o','LineWidth',1.5,'MarkerSize',5)
title('Standard Deviation of Estimation Error','FontSize', 24)
xlabel('Time Steps','FontSize', 24);ylabel('Standard Deviation of Estimation Error','FontSize', 24);
legend('Population','Food Supply','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(7:k,sqrt(data(7:end-1,7)),'--*',7:k,sqrt(data(7:end-1,8)),'-.o')

%% Plot Kalman Estimator Gain
figure
plot(1:k,sqrt(data(2:end,11)),'--*',1:k,sqrt(data(2:end,12)),'-.o','LineWidth',1.5,'MarkerSize',5)
title('Kalman Gain','FontSize', 24)
xlabel('Time Steps','FontSize', 24);ylabel('Kalman Gain','FontSize', 24);
legend('K_1','K_2','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(7:k,sqrt(data(8:end,11)),'--*',7:k,sqrt(data(8:end,12)),'-.o')

%% Find simulation estimation error

err_simulation=[ std(data(:,1)-data(:,3));std(data(:,2)-data(:,4))]; 
fprintf('Standard deviation of estimated error for population based on Pk_plus for %d : timestep is : %d \n' ,k, sqrt(data(end,7)))
fprintf('Standard deviation of estimated error for population from simulation for %d : timestep is : %d \n' ,k, err_simulation(1))
fprintf('Standard deviation of estimated error for food supply based on Pk_plus for %d : timestep is : %d \n' ,k, sqrt(data(end,8)))
fprintf('Standard deviation of estimated error for food supply from simulation for %d : timestep is : %d \n' ,k, err_simulation(2))

%% Plot a priori and a posteriori  for covariance matrix
ti=[0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10];
P_11(1:2:21)=data(1:11,7);P_11(2:2:21)=data(2:11,9);
figure
plot(ti,P_11,0:10,data(:,7),'*',1:10,data(2:11,9),'o','LineWidth',1.5,'MarkerSize',10)
title('Covariance Matrix - a priori and a posteriori population estimation error variance','FontSize', 24)
xlabel('Time Steps','FontSize', 24);ylabel('P_{11}','FontSize', 24);
legend('P_{11} ','P_{11}^+','P_{11}^-','FontSize', 24)

ti=[0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10];
P_22(1:2:21)=data(1:11,8);P_22(2:2:21)=data(2:11,10);
figure
plot(ti,P_22,0:10,data(:,8),'*',1:10,data(2:11,10),'o','LineWidth',1.5,'MarkerSize',10)
title('Covariance Matrix - a priori and a posteriori food supply estimation error variance','FontSize', 24)
xlabel('Time Steps','FontSize', 24);ylabel('P_{22}','FontSize', 24);
legend('P_{22}','P_{22}^+','P_{22}^-' ,'FontSize', 24)
