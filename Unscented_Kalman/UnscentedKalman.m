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
xhat_plus=xhat_plus;                                            % a priori estimate for state x
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

%% Unscented Kalman filter

% Dynamics of system
x0=[300000 -20000 0.001]';  % Initial states                             
rho_0=2;                    % Rho value
g=32.2;                     % Gravitational constant
k=20000;                    % air density altitude relation
Q=0;                        % Process noise
R_k=10000;                  % Measurement noise
M=100000;
a=100000;

xhat_plus=x0;     % a posteriori estimate for state x
P_plus=[1000000 0 0;0 4000000 0;0 0 10]; % a posteriori estimate for covariance P
data=[ x_hat_plus'  diag(P_plus)' 0 0 0 0 0 0]; % Save data
c=2;                                        % Initialize variable
xtilda=chol(3*P_plus);                      % Finding square root  matrix
xhat_tilda=xhat_plus.*ones(3,6);            % xhat vriable to add sigma points
xhat_tilda=xhat_tilda+[xtilda' -xtilda'];   % Sigma points
for i=0.001:0.001:30                        % Loop over time
    xhattilda_dot=[xhat_tilda(2,:); rho_0.*exp(-xhat_tilda(1,:)./k).*xhat_tilda(2,:).^2.*xhat_tilda(3,:)./2-g*ones(1,6);0*ones(1,6)];% Estimate state update
    xhat_tilda=xhat_tilda+xhattilda_dot*0.001;  % Update estimate state
    if mod(i,0.5)==0                            % Checking for time 
        xhat_minus=mean(xhat_tilda,2);          % a priori state estimate
        err=(xhat_tilda-xhat_minus.*ones(3,6)); % Find estimate errors
        P_minus=(err(:,1)*err(:,1)'+err(:,2)*err(:,2)'+err(:,3)*err(:,3)'+err(:,4)*err(:,4)'+err(:,5)*err(:,5)'+err(:,6)*err(:,6)')/6;% a priori error covariance
        yhat_tilda=sqrt(M^2+(xhat_tilda(1,:)-a).^2)';% find output for sigma points   
        yhat=mean(yhat_tilda);                      % Predicted output
        P_y=mean((yhat_tilda-yhat).^2)+R_k;         % Predicted measurement covariance
        err_y=yhat_tilda-yhat;                      % Prediction error
        P_xy=(err(:,1)*err_y(1)+err(:,2)*err_y(2)+err(:,3)*err_y(3)+err(:,4)*err_y(4)+err(:,5)*err_y(5)+err(:,6)*err_y(6))/6;% Cross covariance prediction
        Kk=P_xy*inv(P_y);                           % kalman Gain
        xhat_plus=xhat_minus+Kk*(y(c,1)-yhat);      % a posteriori state estimate
        P_plus=P_minus-Kk*P_y*Kk';                  % a posteriori covariance
        data=[data; xhat_plus' diag(P_plus)' diag(P_minus)' K_k'];     % Saving the data
        xtilda=chol(3*P_plus);                      % Find square root matrix
        xhat_tilda=xhat_plus.*ones(3,6);            % xhat vriable to add sigma points
        xhat_tilda=xhat_tilda+[xtilda' -xtilda'];   % Sigma points
        c=c+1;                          % Incrementing the variable
    end                             
end

% Plot Altitude and Velocity
figure
subplot(2,1,1)
plot(t_m,x_m(:,1), 'LineWidth',2.5);set(gca,'FontSize',20)
title('Altitude','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Altitude','FontSize', 24);
subplot(2,1,2)
plot(t_m,x_m(:,2), 'LineWidth',2.5);set(gca,'FontSize',20)
title('Velocity','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Velocity','FontSize', 24);

% Plot Altitude error
figure
plot(t_m,abs(x_m(:,1)-data_EKF(:,1)),'--*',t_m,abs(x_m(:,1)-data(:,1)),'-.o', 'LineWidth',2.5,'MarkerSize',5);set(gca,'FontSize',20);set(gca, 'YScale', 'log');
title('Altitude Error','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Altitude','FontSize', 24);
legend('Kalman Filter','Unscented Filter','FontSize', 24)

% Plot Velocity error
figure
plot(t_m,abs(x_m(:,2)-data_EKF(:,2)),'--*',t_m,abs(x_m(:,2)-data(:,2)),'-.o', 'LineWidth',2.5,'MarkerSize',5);set(gca,'FontSize',20);set(gca, 'YScale', 'log');
title('Velocity Error','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Velocity','FontSize', 24);
legend('Kalman Filter','Unscented Filter','FontSize', 24)


% Plot Ballistic Coefficient error
figure
plot(t_m,abs(x_m(:,3)-data_EKF(:,3)),'--*',t_m,abs(x_m(:,3)-data(:,3)),'-.o', 'LineWidth',2.5,'MarkerSize',5);set(gca,'FontSize',20);set(gca, 'YScale', 'log');
title('Ballistic Coefficient Error','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Ball. Coeff.','FontSize', 24);
legend('Kalman Filter','Unscented Filter','FontSize', 24)

%% Plot a priori and a posteriori  for covariance matrix for Kalmna filter
ti(3:2:122)=1:60;ti(2:2:121)=1:60;
P_11(1:2:121)=data_EKF(:,7);P_11(2:2:121)=data_EKF(2:end,10);
figure
plot(ti,P_11,0:60,data_EKF(:,7),'*',1:60,data_EKF(2:end,10),'o','LineWidth',1.5,'MarkerSize',5)
title('Covariance Matrix - a priori and a posteriori altitude estimation error variance(KALMAN)','FontSize', 20)
xlabel('Time steps','FontSize', 24);ylabel('P_{11}','FontSize', 24);
legend('P_{11} ','P_{11}^+','P_{11}^-','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(ti(floor(3*end/4):end),P_11(floor(3*end/4):end),45:60,data_EKF(3*end/4:end,7),'*',45:60,data_EKF(3*end/4:end,10),'o','LineWidth',1.5,'MarkerSize',5)

ti(3:2:122)=1:60;ti(2:2:121)=1:60;
P_22(1:2:121)=data_EKF(:,8);P_22(2:2:121)=data_EKF(2:end,11);
figure
plot(ti,P_22,0:60,data_EKF(:,8),'*',1:60,data_EKF(2:end,11),'o','LineWidth',1.5,'MarkerSize',5)
title('Covariance Matrix - a priori and a posteriori velocity estimation error variance(KALMAN)','FontSize', 20)
xlabel('Time steps','FontSize', 24);ylabel('P_{22}','FontSize', 24);
legend('P_{22} ','P_{22}^+','P_{22}^-','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(ti(floor(3*end/4):end),P_22(floor(3*end/4):end),45:60,data_EKF(3*end/4:end,8),'*',45:60,data_EKF(3*end/4:end,11),'o','LineWidth',1.5,'MarkerSize',5)


ti(3:2:122)=1:60;ti(2:2:121)=1:60;
P_33(1:2:121)=data_EKF(:,9);P_33(2:2:121)=data_EKF(2:end,12);
figure
plot(ti,P_33,0:60,data_EKF(:,9),'*',1:60,data_EKF(2:end,12),'o','LineWidth',1.5,'MarkerSize',5)
title('Covariance Matrix - a priori and a posteriori 1/Ballistic coefficient estimation error variance(KALMAN)','FontSize', 20)
xlabel('Time steps','FontSize', 24);ylabel('P_{33}','FontSize', 24);
legend('P_{33} ','P_{33}^+','P_{33}^-','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(ti(floor(3*end/4):end),P_33(floor(3*end/4):end),45:60,data_EKF(3*end/4:end,9),'*',45:60,data_EKF(3*end/4:end,12),'o','LineWidth',1.5,'MarkerSize',5)

%% Plot a priori and a posteriori  for covariance matrix for Kalmna filter
ti(3:2:122)=1:60;ti(2:2:121)=1:60;
P_11(1:2:121)=data(:,4);P_11(2:2:121)=data(2:end,7);
figure
plot(ti,P_11,0:60,data(:,4),'*',1:60,data(2:end,7),'o','LineWidth',1.5,'MarkerSize',5)
title('Covariance Matrix - a priori and a posteriori altitude estimation error variance(UNSCENTED)','FontSize', 20)
xlabel('Time steps','FontSize', 24);ylabel('P_{11}','FontSize', 24);
legend('P_{11} ','P_{11}^+','P_{11}^-','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(ti(floor(3*end/4):end),P_11(floor(3*end/4):end),45:60,data(3*end/4:end,4),'*',45:60,data(3*end/4:end,7),'o','LineWidth',1.5,'MarkerSize',5)

ti(3:2:122)=1:60;ti(2:2:121)=1:60;
P_22(1:2:121)=data(:,5);P_22(2:2:121)=data(2:end,8);
figure
plot(ti,P_22,0:60,data(:,5),'*',1:60,data(2:end,8),'o','LineWidth',1.5,'MarkerSize',5)
title('Covariance Matrix - a priori and a posteriori velocity estimation error variance(UNSCENTED)','FontSize', 20)
xlabel('Time steps','FontSize', 24);ylabel('P_{22}','FontSize', 24);
legend('P_{22} ','P_{22}^+','P_{22}^-','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(ti(floor(3*end/4):end),P_22(floor(3*end/4):end),45:60,data(3*end/4:end,5),'*',45:60,data(3*end/4:end,8),'o','LineWidth',1.5,'MarkerSize',5)


ti(3:2:122)=1:60;ti(2:2:121)=1:60;
P_33(1:2:121)=data(:,6);P_33(2:2:121)=data(2:end,9);
figure
plot(ti,P_33,0:60,data(:,6),'*',1:60,data(2:end,9),'o','LineWidth',1.5,'MarkerSize',5)
title('Covariance Matrix - a priori and a posteriori 1/Ballistic coefficient estimation error variance(UNSCENTED)','FontSize', 20)
xlabel('Time steps','FontSize', 24);ylabel('P_{33}','FontSize', 24);
legend('P_{33} ','P_{33}^+','P_{33}^-','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(ti(floor(3*end/4):end),P_33(floor(3*end/4):end),45:60,data(3*end/4:end,6),'*',45:60,data(3*end/4:end,9),'o','LineWidth',1.5,'MarkerSize',5)


