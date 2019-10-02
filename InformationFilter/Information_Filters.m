% Program by Sreevalsan S Menon(sm2hm@mst.edu)

clc
clear 
close all
rng default
%% System Paramters
t = 0:0.001:5;  % Simulation Time
F = 0.001;      % Coefficient of viscous friction
J = 0.00018;    % Moment of inertia
I = 0.003;      % Winding inductance
lambda = 0.1;   % Motor constant
Re = 1.9;       % Winding resistance
x = zeros(4,1); % Initial states
%% Noise matrix
R = [10e-6 0; 0 10e-6];                                 % Measurement noise
Q = [11.11 0 0 0;0 11.11 0 0;0 0 0.25 0;0 0 0 1e-6];    % Process Noise


%% Extended Information Filter  
Y=1/1000*eye(4);                                        % Information matrix initialization
x_hat=[0.1 0.2 0.3 0.4]';                               % Initial state estimate

%% Cubature Information Filter
Pc=inv(1000*eye(4));                                    % Covariance Iniatilzation
Sc=chol(Pc,'lower');                                    % Finding square root of covriance matrix
Syc=inv(Sc);                                            % Finding square root information matrix
x_hatc=[0.1 0.2 0.3 0.4]';                              % Initial state estimate
Xi =  repmat(x_hatc,1,8) + Sc*sqrt(4)*[eye(4) -eye(4)]; % Finding cubature points

%% Square root Kalman Information Filter
Psk=1000*eye(4);                                        % Covariance Iniatilzation
Ysk=chol(Psk\eye(4));                                   % Finding square root information matrix
x_hatsk=[0.1 0.2 0.3 0.5]';                             % Initial state estimate

%% Square root Cubature Information Filter
Psc=1000*eye(4);                                        % Covariance Iniatilzation
Qsqrt=chol(Q);                                          % Square root of process noise
Ssc=chol(Psc\eye(4),'lower');                           % Square root of covariance matrix
Sysc=inv(Ssc);                                          % Square root information matrix
x_hatsc=[0.1 0.2 0.3 0.5]';                             % Initial state estimate
Xis =  repmat(x_hatsc,1,8) + Ssc*sqrt(4)*[eye(4) -eye(4)]; % Finding cubature points
%% Data array
data=[];              % Array to save loop values

%% Prediction and Measurement Update
for i=1:size(t,2)
    %% Actual System
    u1=sin(0.002*pi*i);                                         % State Input u1
    u2=cos(0.002*pi*i);                                         % State Input u2
    x=[ x(1)+0.001*(-Re*x(1)/I + x(3)*lambda*sin(x(4))/I+u1/I)
        x(2)+0.001*(-Re*x(2)/I - x(3)*lambda*cos(x(4))/I+u2/I)
        x(3)+0.001*(-3*lambda*x(1)*sin(x(4))/(2*J)+3*lambda*x(2)*cos(x(4))/(2*J)-F*x(3)/J)
        x(4)+0.001*x(3)];                                       % State update
    x=x+0.001*sqrt(Q)*randn(4,1);                               % Add Process noise
    y=x(1:2)+sqrt(R)*randn(2,1);                                % Finding output 
    
    %% Extended Information Filter(EIF)
    % Prediction
    x_hat=[ x_hat(1)+0.001*(-Re*x_hat(1)/I + x_hat(3)*lambda*sin(x_hat(4))/I+u1/I)
            x_hat(2)+0.001*(-Re*x_hat(2)/I - x_hat(3)*lambda*cos(x_hat(4))/I+u2/I)
            x_hat(3)+0.001*(-3*lambda*x_hat(1)*sin(x_hat(4))/(2*J)+3*lambda*x_hat(2)*cos(x_hat(4))/(2*J)-F*x_hat(3)/J)
            x_hat(4)+0.001*x_hat(3)];                           % State estimate update EIF
    A=[-Re/I 0 lambda/I*sin(x_hat(4)) x_hat(3)*lambda*cos(x_hat(4))/I
        0 -Re/I -lambda*cos(x_hat(4))/I x_hat(3)*lambda*sin(x_hat(4))/I
       -3/2*lambda*sin(x_hat(4))/J 3*lambda*cos(x_hat(4))/(2*J) -F/J -3*lambda*(x_hat(1)*cos(x_hat(4))/(2*J)+x_hat(2)*sin(x_hat(4)))
        0 0 1 0];                                               % Linearization of stae matrix
    Y=inv(A*(Y\A')+Q);                                          % Update information matrix
    y_hat=Y*x_hat;                                              % Update information states
     
    % Measurement Update
    C=[1 0 0 0;0 1 0 0];                                        % Output matrix
    v_k=y-[x_hat(1);x_hat(2)];                                  % Measurement residual
    i_k=C'*(R\(v_k+C*x_hat));                                   % Information state contribution
    I_k=C'*(R\C);                                               % Information matrix contrbution
    y_hat=y_hat+i_k;                                            % Update Information state
    Y=Y+I_k;                                                    % Update Information matrix
    x_hat=Y\y_hat;                                              % Update state estimate
    
    %% Cubature information Filter
    %Prediction
    Xi=[ Xi(1,:)+0.001.*(-Re.*Xi(1,:)./I + Xi(3,:).*lambda.*sin(Xi(4,:))./I+u1/I*ones(1,8))
         Xi(2,:)+0.001.*(-Re*Xi(2,:)./I - Xi(3,:).*lambda.*cos(Xi(4,:))./I+u2/I*ones(1,8))
         Xi(3,:)+0.001.*(-3.*lambda.*Xi(1,:).*sin(Xi(4,:))./(2.*J)+3.*lambda.*Xi(2,:).*cos(Xi(4,:))./(2.*J)-F.*Xi(3,:)./J)
         Xi(4,:)+0.001*Xi(3,:)];                                % Update cubature estimates
    x_hatc=mean(Xi,2);                                          % Find mean estimate
    Pc=(1/8)*(Xi*Xi')-x_hatc*x_hatc'+Q;                         % Predicted covariance matrix
    Sc=chol(Pc);                                                % Square root of covariance matrix
    Yc=inv(Pc);                                                 % Predicted information matrix
    yhat_c=Yc*x_hatc;                                           % Predicted information states
    
    % Measurement Update
    Xi =  repmat(x_hatc,1,8) + Sc*sqrt(4)*[eye(4) -eye(4)];     % Update cubature points
    z_i=Xi(1:2,:);                                              % Output estimate
    z_k=mean(z_i,2);                                            % Output estimate mean
    vk=y-z_k;                                                   % Measurement residual
    Pxz=(1/8)*Xi*z_i'-x_hatc*z_k';                              % Error cross covariance matrix
    I_kc=Yc*Pxz*inv(R)*Pxz'*Yc';                                % Information matrix contrbution
    i_kc=Yc*Pxz*inv(R)*(vk+Pxz'*Yc'*x_hatc);                    % Information state contrbution
    Yc=Yc+I_kc;                                                 % Update Information matrix
    yhat_c=yhat_c+i_kc;                                         % Update Information states
    x_hatc=Yc\yhat_c;                                           % Update state estimate
    Pc=Yc\eye(4);                                               % Update covariance matrix
    Sc=chol(diag(diag(Pc)),'lower');                            % Find square root of covariance matrix
    Xi =  repmat(x_hatc,1,8) + Sc*sqrt(4)*[eye(4) -eye(4)];     % Update cubate points
    
% %% Square root extended information filter
%     % Prediction
%     x_hatsk=[ x_hatsk(1)+0.001*(-Re*x_hatsk(1)/I + x_hatsk(3)*lambda*sin(x_hatsk(4))/I+u1/I)
%             x_hatsk(2)+0.001*(-Re*x_hatsk(2)/I - x_hatsk(3)*lambda*cos(x_hatsk(4))/I+u2/I)
%             x_hatsk(3)+0.001*(-3*lambda*x_hatsk(1)*sin(x_hatsk(4))/(2*J)+3*lambda*x_hatsk(2)*cos(x_hatsk(4))/(2*J)-F*x_hatsk(3)/J)
%             x_hatsk(4)+0.001*x_hatsk(3)];                       % Predict state estimate
%     A=[-Re/I 0 lambda/I*sin(x_hatsk(4)) x_hatsk(3)*lambda*cos(x_hatsk(4))/I
%         0 -Re/I -lambda*cos(x_hatsk(4))/I x_hatsk(3)*lambda*sin(x_hatsk(4))/I
%        -3/2*lambda*sin(x_hatsk(4))/J 3*lambda*cos(x_hatsk(4))/(2*J) -F/J -3*lambda*(x_hatsk(1)*cos(x_hatsk(4))/(2*J)+x_hatsk(2)*sin(x_hatsk(4)))
%         0 0 1 0];                                               % Linearization of stae matrix
%     Ysk=chol((inv(A*((Ysk*Ysk')\A')+Q)),'lower');               % Predict square root of information matrix
%     y_hatsk=Ysk*x_hatsk;                                        % Predict information states
%     
%     % Measurement Update
%     C=[0 1 0 0;1 0 0 0];                                        % Output matrix
%     [W, T]  = qr([Ysk C'*sqrt(R);y_hatsk' y'*sqrt(R)]',0);T=T';  % Householder transformation QR decomposition
%     Ysk=tril(T(1:4,1:4));                                       % Update information matrix
%     y_hatsk=T(5,1:4)';                                          % Update information state
%     x_hatsk=Ysk\y_hatsk;                                        % Update state estimate
%     
%% Square root cubature information filter
    %Prediction
    Xis=[ Xis(1,:)+0.001.*(-Re.*Xis(1,:)./I + Xis(3,:).*lambda.*sin(Xis(4,:))./I+u1/I*ones(1,8))
         Xis(2,:)+0.001.*(-Re*Xis(2,:)./I - Xis(3,:).*lambda.*cos(Xis(4,:))./I+u2/I*ones(1,8))
         Xis(3,:)+0.001.*(-3.*lambda.*Xis(1,:).*sin(Xis(4,:))./(2.*J)+3.*lambda.*Xis(2,:).*cos(Xis(4,:))./(2.*J)-F.*Xis(3,:)./J)
         Xis(4,:)+0.001*Xis(3,:)];                                  % Update cubature estimates    
    x_hatsc=mean(Xis,2);                                            % Find mean estimate
    Xt=(Xis-repmat(x_hatsc,1,8))/sqrt(8);                           % Weighted center matrix
    [~,Ssc]=qr([Xt Qsqrt]',0);Ssc=Ssc';                             % Predict square root of covariance
    Sysc=eye(4)/Ssc';                                               % Predict square root of information matrix
    yhat_sc=Sysc*Sysc'*x_hatsc;                                     % Predict information states
    
    % Measurement Update
    Xis =  repmat(x_hatsc,1,8) + Ssc*sqrt(4)*[eye(4) -eye(4)];      % Update cubature points
    z_is=Xis(1:2,:);                                                % Output estimate
    z_ks=mean(z_is,2);                                              % Mean output estimate
    Xs=(Xis-repmat(x_hatsc,1,8))/sqrt(8);                           % Weighted center matrix for X
    Z=(z_is-repmat(z_ks,1,8))/sqrt(8);                              % Weighted center matrix for Z
    [~,Ss]=qr([Z sqrt(R); Xs zeros(4,2)]',0);Ss=Ss';                % Perform QR decomposition
    Xs=Ss(1:2,1:2);Ys=Ss(3:end,1:2);                                % 
    I_sc=Sysc*Sysc'*(Ys*Xs')*sqrt(pinv(R));                         % Information matrix contribution
    i_sc=(I_sc*sqrt(pinv(R)))*(y-z_ks)+I_sc*I_sc'*x_hatsc;          % Information state contribution
    [~,Sysc]=qr([Sysc I_sc]',0);Sysc=Sysc';                         % Update square root of Information matrix
    yhat_sc=yhat_sc+i_sc;                                           % Update Information state
    Ssc=eye(4)/Sysc';                                               % Update square root of covariance matrix
    x_hatsc=Ssc*Ssc'*yhat_sc;                                       % Find state estimate
    Xis =  repmat(x_hatsc,1,8) + Ssc*sqrt(4)*[eye(4) -eye(4)];      % Update cubature points
    
    %% Saving estimates
    if mod(t(i),0.01)==0
       data=[data; x' x_hat' x_hatc' x_hatsk' x_hatsc'];        % Saving the data
    end  
end
%% Plot X1

figure
plot(0:0.01:5,data(:,1),0:0.01:5,data(:,5),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_1 in Amps','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_1 in Amps','FontSize', 24);
legend('ACTUAL','EIF','FontSize', 24)
axes('Position',[.65 .2 .2 .2])
box on
plot(3.8:0.01:4,data(381:401,1),3.8:0.01:4,data(381:401,5),'--*', 'LineWidth',1.5,'MarkerSize',2.5)

figure
plot(0:0.01:5,data(:,1),0:0.01:5,data(:,9),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_1 in Amps','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_1 in Amps','FontSize', 24);
legend('ACTUAL','CIF','FontSize', 24)
axes('Position',[.65 .2 .2 .2])
box on
plot(3.8:0.01:4,data(381:401,1),3.8:0.01:4,data(381:401,9),'--*', 'LineWidth',1.5,'MarkerSize',2.5)


figure
plot(0:0.01:5,data(:,1),0:0.01:5,data(:,17),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_1 in Amps','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_1 in Amps','FontSize', 24);
legend('ACTUAL','SCIF','FontSize', 24)
axes('Position',[.65 .2 .2 .2])
box on
plot(3.8:0.01:4,data(381:401,1),3.8:0.01:4,data(381:401,17),'--*', 'LineWidth',1.5,'MarkerSize',2.5)

%% Plot X2

figure
plot(0:0.01:5,data(:,2),0:0.01:5,data(:,6),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_2 in Amps','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_2 in Amps','FontSize', 24);
legend('ACTUAL','EIF','FontSize', 24)
axes('Position',[.65 .2 .2 .2])
box on
plot(3.8:0.01:4,data(381:401,2),3.8:0.01:4,data(381:401,6),'--*', 'LineWidth',1.5,'MarkerSize',2.5)

figure
plot(0:0.01:5,data(:,2),0:0.01:5,data(:,10),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_2 in Amps','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_2 in Amps','FontSize', 24);
legend('ACTUAL','CIF','FontSize', 24)
axes('Position',[.65 .2 .2 .2])
box on
plot(3.8:0.01:4,data(381:401,2),3.8:0.01:4,data(381:401,10),'--*', 'LineWidth',1.5,'MarkerSize',2.5)


figure
plot(0:0.01:5,data(:,2),0:0.01:5,data(:,18),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_2 in Amps','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_2 in Amps','FontSize', 24);
legend('ACTUAL','SCIF','FontSize', 24)
axes('Position',[.65 .2 .2 .2])
box on
plot(3.8:0.01:4,data(381:401,2),3.8:0.01:4,data(381:401,18),'--*', 'LineWidth',1.5,'MarkerSize',2.5)


%% Plot X3
figure
plot(0:0.01:5,data(:,3),0:0.01:5,data(:,7),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_3 in Rad/s','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_3 in Rad/s','FontSize', 24);
legend('ACTUAL','EIF','FontSize', 24)
axes('Position',[.65 .2 .2 .2])
box on
plot(3.8:0.01:4,data(381:401,3),3.8:0.01:4,data(381:401,7),'--*', 'LineWidth',1.5,'MarkerSize',2.5)

figure
plot(0:0.01:5,data(:,3),0:0.01:5,data(:,11),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_3 in Rad/s','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_3 in Rad/s','FontSize', 24);
legend('ACTUAL','CIF','FontSize', 24)
axes('Position',[.65 .6 .2 .2])
box on
plot(3.8:0.01:4,data(381:401,3),3.8:0.01:4,data(381:401,11),'--*', 'LineWidth',1.5,'MarkerSize',2.5)


figure
plot(0:0.01:5,data(:,3),0:0.01:5,data(:,19),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_3 in Rad/s','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_3 in Rad/s','FontSize', 24);
legend('ACTUAL','SCIF','FontSize', 24)
axes('Position',[.65 .6 .2 .2])
box on
plot(3.8:0.01:4,data(381:401,3),3.8:0.01:4,data(381:401,19),'--*', 'LineWidth',1.5,'MarkerSize',2.5)

%% Plot X4
figure
plot(0:0.01:5,mod(data(:,4),2*pi),0:0.01:5,mod(data(:,8),2*pi),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_4 in Rad','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_4 in Rad','FontSize', 24);
legend('ACTUAL','EIF','FontSize', 24)
axes('Position',[.65 .2 .2 .2])
box on
plot(3.8:0.01:4,mod(data(381:401,4),2*pi),3.8:0.01:4,mod(data(381:401,8),2*pi),'--*', 'LineWidth',1.5,'MarkerSize',2.5)

figure
plot(0:0.01:5,mod(data(:,4),2*pi),0:0.01:5,mod(data(:,12),2*pi),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_4 in Rad','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_4 in Rad','FontSize', 24);
legend('ACTUAL','CIF','FontSize', 24)
axes('Position',[.65 .2 .2 .2])
box on
plot(3.8:0.01:4,mod(data(381:401,4),2*pi),3.8:0.01:4,mod(data(381:401,12),2*pi),'--*', 'LineWidth',1.5,'MarkerSize',2.5)


figure
plot(0:0.01:5,mod(data(:,4),2*pi),0:0.01:5,mod(data(:,20),2*pi),'--*', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('X_4 in Rad','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('X_4 in Rad','FontSize', 24);
legend('ACTUAL','SCIF','FontSize', 24)
axes('Position',[.65 .2 .2 .2])
box on
plot(3.8:0.01:4,mod(data(381:401,4),2*pi),3.8:0.01:4,mod(data(381:401,20),2*pi),'--*', 'LineWidth',1.5,'MarkerSize',2.5)

%% Plot RMSE
figure
plot(0:0.01:5,abs(data(:,1)-data(:,5)),0:0.01:5,abs(data(:,1)-data(:,9)), 0:0.01:5,abs(data(:,1)-data(:,17)), 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('Error X1','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Error X1','FontSize', 24);
legend('EIF','CIF','SCIF','FontSize', 24)

figure
plot(0:0.01:5,abs(data(:,2)-data(:,6)),0:0.01:5,abs(data(:,2)-data(:,10)), 0:0.01:5,abs(data(:,2)-data(:,18)), 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('Error X2','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Error X2','FontSize', 24);
legend('EIF','CIF','SCIF','FontSize', 24)

figure
plot(0:0.01:5,abs(data(:,3)-data(:,7)),'m',0:0.01:5,abs(data(:,3)-data(:,11)),'r', 0:0.01:5,abs(data(:,3)-data(:,19)),'k', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('Error X3','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Error X3','FontSize', 24);
legend('EIF','CIF','SCIF','FontSize', 24)
axes('Position',[.65 .6 .2 .2])
box on
plot(3.8:0.01:5,abs(data(381:end,3)-data(381:end,11)),'r' , 'LineWidth',1.5,'MarkerSize',2.5)
axes('Position',[.65 0.3 .2 .2])
box on
plot(3.8:0.01:5,abs(data(381:end,3)-data(381:end,19)), 'k','LineWidth',1.5,'MarkerSize',2.5)

figure
plot(0:0.01:5,abs(data(:,4)-data(:,8)),'m',0:0.01:5,abs(data(:,4)-data(:,12)),'r', 0:0.01:5,abs(data(:,4)-data(:,20)),'k', 'LineWidth',1.5,'MarkerSize',2.5); 
set(gca,'FontSize',20);
title('Error X4','FontSize', 24)
xlabel('Time ','FontSize', 24);ylabel('Error X4','FontSize', 24);
legend('EIF','CIF','SCIF','FontSize', 24)
axes('Position',[.65 .6 .2 .2])
box on
plot(3.8:0.01:5,abs(data(381:end,4)-data(381:end,12)),'r' , 'LineWidth',1.5,'MarkerSize',2.5)
axes('Position',[.65 0.3 .2 .2])
box on
plot(3.8:0.01:5,abs(data(381:end,4)-data(381:end,20)), 'k','LineWidth',1.5,'MarkerSize',2.5)