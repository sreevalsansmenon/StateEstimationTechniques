% Program by Sreevalsan S Menon(sm2hm@mst.edu)
clc                                                                         % Clear command window
clear                                                                       % clear workspace
% rng('default')                                                              % For reproducibility
close all                                                                   % Close all figures

i=1;                                    % Initializing i
T=0.1;                                  % Step time
% System dynamics 
Q=diag([0,0,4,4]);                      % Process covariance matrix
R=diag([1,1]);                          % Measurement covariance matrix
N1=20;E1=0;N2=0;E2=20;                  % Tracking station location
x=[0;0;50;50];                          % Initial states
F=[1 0 T 0;0 1 0 T;0 0 1 0;0 0 0 1];    % System matrix
L=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];    % Process covariance 

xhat_plus=x;                            % Initializing a posteriori state estimate
P_plus=0*eye(4);                        % initializing a posteriori covariance matrix P
data=[x' xhat_plus' xhat_plus' diag(P_plus)' 0 0 0 0]; % Saving data

% Finding the estimate
for t=0:T:60                            % Iterating over time
    
    x=F*x+[0 0 2*randn 2*randn]';       % Finding state
    y(:,i)=[sqrt((x(1)-N1)^2+((x(2)-E1)^2))+randn;sqrt((x(1)-N2)^2+((x(2)-E2)^2))+randn]; % Finding system output
    
    P_minus=F*P_plus*F'+L*Q*L';         % Updating a priori covariance matrix P
    xhat_minus=F*xhat_plus;             % Updating a priori state matrix estimate
    
    % Updating output matrix at nominal values
    H=[(xhat_minus(1)-N1)/sqrt((xhat_minus(1)-N1)^2+((xhat_minus(2)-E1)^2)) (xhat_minus(2)-E1)/sqrt((xhat_minus(1)-N1)^2+((xhat_minus(2)-E1)^2)) 0 0;
       (xhat_minus(1)-N2)/sqrt((xhat_minus(1)-N2)^2+((xhat_minus(2)-E2)^2)) (xhat_minus(2)-E2)/sqrt((xhat_minus(1)-N2)^2+((xhat_minus(2)-E2)^2)) 0 0];
    
    M=eye(2);                           % Measurement covariance matrix 
    
    K=P_minus*H'*inv(H*P_minus*H'+M*R*M); % Kalman estmate gain
    xhat_plus=xhat_minus+K*(y(:,i)-[sqrt((xhat_minus(1)-N1)^2+(xhat_minus(2)-E1)^2);sqrt((xhat_minus(1)-N2)^2+(xhat_minus(2)-E2)^2)]);% Updating a posteriori covariance matrix P
    P_plus=(eye(4)-K*H)*P_minus;        % Updating a posteriori state matrix estimate
    i=i+1;                              % incrementing i variable
    data=[data; x' xhat_plus' xhat_minus' diag(P_plus)' diag(P_minus)']; % Saving data
end

% Plot North coordinate Position error
figure
plot(0:T:60,abs(data(2:end,1)-data(2:end,5)),'LineWidth',1.5);set(gca,'FontSize',20);
title('North coordinate Position Error','FontSize', 24)
xlabel('Time','FontSize', 24);ylabel('Estimation Error','FontSize', 24);

% Plot East coordinate Position error
figure
plot(0:T:60,abs(data(2:end,2)-data(2:end,6)),'LineWidth',1.5);set(gca,'FontSize',20);
title('East coordinate Position Error','FontSize', 24)
xlabel('Time','FontSize', 24);ylabel('Estimation Error','FontSize', 24);

% North coordinate Velocity Error
figure
plot(0:T:60,abs(data(2:end,3)-data(2:end,7)),'LineWidth',1.5);set(gca,'FontSize',20);
title('North coordinate Velocity Error','FontSize', 24)
xlabel('Time','FontSize', 24);ylabel('Estimation Error','FontSize', 24);

% East coordinate Velocity Error
figure
plot(0:T:60,abs(data(2:end,4)-data(2:end,8)),'LineWidth',1.5);set(gca,'FontSize',20);
title('East coordinate Velocity Error','FontSize', 24)
xlabel('Time','FontSize', 24);ylabel('Estimation Error','FontSize', 24);

%% Find simulation estimation error

err_simulation=[ std(data(2:end,1)-data(2:end,5));std(data(2:end,2)-data(2:end,6));std(data(2:end,3)-data(2:end,7));std(data(2:end,4)-data(2:end,8))]; 
fprintf('Standard deviation of estimated error for North coordinate Position from simulation is %d : and theoretical is : %d \n' ,err_simulation(1),sqrt(P_plus(1,1)))
fprintf('Standard deviation of estimated error for East coordinate Position from simulation is %d : and theoretical is : %d \n' ,err_simulation(2),sqrt(P_plus(2,2)))
fprintf('Standard deviation of estimated error for North coordinate Velocity from simulation is %d : and theoretical is : %d \n' ,err_simulation(3),sqrt(P_plus(3,3)))
fprintf('Standard deviation of estimated error for East coordinate Velocity from simulation is %d : and theoretical is : %d \n' ,err_simulation(4),sqrt(P_plus(4,4)))


%% Plot a priori and a posteriori  for covariance matrix
ti(3:2:602*2)=1:601;ti(2:2:602*2-1)=1:601;
P_11(1:2:1203)=data(:,13);P_11(2:2:1203)=data(2:end,17);
figure
plot(ti(end-50:end),P_11(end-50:end),576:601,data(577:602,13),'*',577:601,data(578:602,17),'o','LineWidth',1.5,'MarkerSize',10)
title('Covariance Matrix - a priori and a posteriori North coordinate Position estimation error variance(Hybrid)','FontSize', 24)
xlabel('Time steps','FontSize', 24);ylabel('P_{11}','FontSize', 24);
legend('P_{11} ','P_{11}^+','P_{11}^-','FontSize', 24)

ti(3:2:602*2)=1:601;ti(2:2:602*2-1)=1:601;
P_22(1:2:1203)=data(:,14);P_22(2:2:1203)=data(2:end,18);
figure
plot(ti(end-50:end),P_22(end-50:end),576:601,data(577:602,14),'*',577:601,data(578:602,18),'o','LineWidth',1.5,'MarkerSize',10)
title('Covariance Matrix - a priori and a posteriori East coordinate Position estimation error variance(Hybrid)','FontSize', 24)
xlabel('Time steps','FontSize', 24);ylabel('P_{22}','FontSize', 24);
legend('P_{22} ','P_{22}^+','P_{22}^-','FontSize', 24)

ti(3:2:602*2)=1:601;ti(2:2:602*2-1)=1:601;
P_33(1:2:1203)=data(:,15);P_33(2:2:1203)=data(2:end,19);
figure
plot(ti(end-50:end),P_33(end-50:end),576:601,data(577:602,15),'*',577:601,data(578:602,19),'o','LineWidth',1.5,'MarkerSize',10)
title('Covariance Matrix - a priori and a posteriori North coordinate Velocity estimation error variance(Hybrid)','FontSize', 24)
xlabel('Time steps','FontSize', 24);ylabel('P_{33}','FontSize', 24);
legend('P_{33} ','P_{33}^+','P_{33}^-','FontSize', 24)

ti(3:2:602*2)=1:601;ti(2:2:602*2-1)=1:601;
P_44(1:2:1203)=data(:,16);P_44(2:2:1203)=data(2:end,20);
figure
plot(ti(end-50:end),P_44(end-50:end),576:601,data(577:602,16),'*',577:601,data(578:602,20),'o','LineWidth',1.5,'MarkerSize',10)
title('Covariance Matrix - a priori and a posteriori East coordinate Velocity estimation error variance(Hybrid)','FontSize', 24)
xlabel('Time steps','FontSize', 24);ylabel('P_{44}','FontSize', 24);
legend('P_{44} ','P_{44}^+','P_{44}^-','FontSize', 24)



