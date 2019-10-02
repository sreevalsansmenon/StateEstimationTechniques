% Program by Sreevalsan S Menon(sm2hm@mst.edu)

clc                         % Clear command window
clear                       % clear workspace
close all                   % Close all figures
rng default


Q_k=10;                      % Process noise
R_k=1;                      % Measurement noise
x_k=0.1;                    % Initial state
xk_hat=0.1;                 % Initial state estimate
Pk_plus=2;                  % Initial Covariance matrix
kpdf=[];                    % Variable to save kalman pdf
data_extended=[ x_k xk_hat Pk_plus 0 0 0];  % Saving data

N=500;                      % Setting number of particles
x_pk=x_k+sqrt(2)*normrnd(0,1,[N,1]);% Generating N particles
ppdf=[];                    % Variable to save particle pdf
data_particle=[xk_hat x_pk' zeros(1,N)];% Save data

for i=1:50                          
        % Actual system
        x_k=0.5*x_k+25*x_k/(1+x_k^2)+8*cos(1.2*(i-1))+sqrt(Q_k)*randn;% State update
        y_k=x_k^2/20+sqrt(R_k)*randn;   % Measurement update
        
        %% Extended Kalman Filter 
        F_k=25*x_k/(x_k^2+1)+x_k/2;     % System matrix
        L_k=1;                          % Process covariance

        Pk_minus=F_k*Pk_plus*F_k'+L_k*Q_k*L_k'; % a priori covariance
        xk_minus=0.5*xk_hat+25*xk_hat/(1+xk_hat^2)+8*cos(1.2*(i-1));% a priori state estimate
        
        H_k=xk_minus/10;                % Partial derivative of measurement function
        M_k=1;                          % Partial derivative of measurement variance
        
        K_k=Pk_minus*H_k'*inv(H_k*Pk_minus*H_k'+M_k*R_k*M_k');% Updating Kalman estimator gain
        xk_hat=xk_hat+K_k*(y_k-xk_hat^2/20);% Finding state estimate
        Pk_plus=(1-K_k*H_k)*Pk_minus;       % a posteriori covariance
        data_extended=[data_extended;x_k xk_hat Pk_plus Pk_minus K_k y_k]; % Saving the data
        
        %% Particle filter
        x_pk_minus=0.5.*x_pk+25.*x_pk./(1+x_pk.^2)+8*cos(1.2*(i-1))+sqrt(Q_k)*normrnd(0,1,[N,1]); % a priori estimate for state
        y_pk=x_pk_minus.^2/20;              % Measurement estimate
        q_i=1/((2*pi)^0.5*R_k^0.5)*exp(-(y_k-y_pk).*inv(R_k).*(y_k-y_pk)/2); % Compute relative likelihood
        q_i=q_i/sum(q_i);                   % Weight relative likelihood
        
        % Resampling based on likelihood
        for k=1:N
            r=rand;                         % Find random number between 0 and 1
            rsum=0;                         % Temporary variable to add cumulative probability value
            for j=1:N
                rsum=rsum+q_i(j);           % Finding cumulative probability
                if rsum>=r                  % Check cumulative probability greater than r
                    xpk(k)=x_pk_minus(j);   % Save the a posteriori
                    break;
                end
            end
        end
        
        x_pk=xpk';                          % a posteriori particles
        xpk_hat=mean(x_pk);                 % algebraic mean of particle
        data_particle=[data_particle;xpk_hat x_pk' q_i']; % Save data

        % Particle filter pdf
        pdf= zeros(41,1);
        for m = -20 : 20
            for l = 1 : N
                if (m<=x_pk(l))&&(x_pk(l)<m+1)
                    pdf(m+21)=pdf(m+21)+1;
                end
            end
        end
        ppdf=[ppdf;pdf'/N];
        % Kalman filter pdf
        m=-60:60;
        pdf = 1/(sqrt(Pk_plus)*sqrt(2*pi)) .* exp(-(m - xk_hat).^2 / 2 / Pk_plus);
        kpdf=[kpdf;pdf];
end

%% Plot Particle Posterior PDF at Fixed Time Points

figure
sgtitle('Posterior PDF at Fixed Time Points(N=500; Dashed red curve = true state value)','FontSize', 24)
subplot(3,3,1)
bar(-20:20,ppdf(5,:));y1=get(gca,'ylim');hold on;plot([data_extended(6,1) data_extended(6,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=5','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,2)
bar(-20:20,ppdf(10,:));y1=get(gca,'ylim');hold on;plot([data_extended(11,1) data_extended(11,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=10','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,3)
bar(-20:20,ppdf(15,:));y1=get(gca,'ylim');hold on;plot([data_extended(16,1) data_extended(16,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=15','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,4)
bar(-20:20,ppdf(20,:));y1=get(gca,'ylim');hold on;plot([data_extended(21,1) data_extended(21,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=20','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,5)
bar(-20:20,ppdf(25,:));y1=get(gca,'ylim');hold on;plot([data_extended(26,1) data_extended(26,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=25','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,6)
bar(-20:20,ppdf(30,:));y1=get(gca,'ylim');hold on;plot([data_extended(31,1) data_extended(31,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=30','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,7)
bar(-20:20,ppdf(35,:));y1=get(gca,'ylim');hold on;plot([data_extended(36,1) data_extended(36,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=35','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,8)
bar(-20:20,ppdf(40,:));y1=get(gca,'ylim');hold on;plot([data_extended(41,1) data_extended(41,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=40','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,9)
bar(-20:20,ppdf(45,:));y1=get(gca,'ylim');hold on;plot([data_extended(46,1) data_extended(46,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=45','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);


%%  Particle Evolution of Posterior PDF
figure
patchColor = [0 0.4470 0.7410];
[X,Y] = meshgrid(1:2:50,-20:20);
hFill = fill3(X,Y, ppdf(1:2:end,:), patchColor, 'LineWidth', 1, 'EdgeColor', patchColor,'FaceAlpha', 0.5);
grid on
hold on
plot(1:2:50,data_extended(2:2:end,1),'r','LineWidth',4)
title('Particle Evolution of Posterior PDF')
xlabel('K','FontSize', 24);ylabel('X','FontSize', 24);zlabel('PDF','FontSize', 24);set(gca,'FontSize',20);ylim([-20 20])

%% True & Estimated State With 5 & 95 Percentiles

figure
sgtitle('Particle Filter -TRUE AND ESTIMATED STATE WITH 5 & 95 PERCENTILE BOUNDS','FontSize', 24)
plot(1:50,data_extended(2:end,1),'b', 'LineWidth',3); hold on 
plot(1:50,data_particle(2:end,1),'--or','LineWidth',3)
plot(1:50,prctile( data_particle(2:end,2:N+1)', 95 ),'--g',1:50,prctile( data_particle(2:end,2:N+1)', 5 ),'--g', 'LineWidth',3)
legend('Truth','Estimate','5/95 Percentiles')
xlabel('K','FontSize', 24);ylabel('X','FontSize', 24);

%% Particle Estimation Error and 2-Sigma Bounds

figure 
sgtitle('Particle Filter -ESTIMATION ERROR AND 2-SIGMA BOUNDS','FontSize', 24)
plot(1:50,data_extended(2:end,1)-data_particle(2:end,1),'b','LineWidth',3);hold on;set(gca,'FontSize',20);
plot(1:50, 2.*std(data_particle(2:end,2:N+1)'),'r',1:50,-2.*std(data_particle(2:end,2:N+1)'),'r','LineWidth',3)
legend('Error','Bound')
xlabel('K','FontSize', 24);ylabel('ERROR','FontSize', 24);

%% Plot Kalman Particle Posterior PDF at Fixed Time Points
figure
sgtitle('EKF Posterior PDF at Fixed Time Points(Dashed red curve = true state value)','FontSize', 24)
subplot(3,3,1)
plot(-60:60,kpdf(5,:),'LineWidth',3);y1=get(gca,'ylim');hold on;plot([data_extended(6,1) data_extended(6,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=5','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,2)
plot(-60:60,kpdf(10,:),'LineWidth',3);y1=get(gca,'ylim');hold on;plot([data_extended(11,1) data_extended(11,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=10','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,3)
plot(-60:60,kpdf(15,:),'LineWidth',3);y1=get(gca,'ylim');hold on;plot([data_extended(16,1) data_extended(16,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=15','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,4)
plot(-60:60,kpdf(20,:),'LineWidth',3);y1=get(gca,'ylim');hold on;plot([data_extended(21,1) data_extended(21,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=20','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,5)
plot(-60:60,kpdf(25,:),'LineWidth',3);y1=get(gca,'ylim');hold on;plot([data_extended(26,1) data_extended(26,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=25','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,6)
plot(-60:60,kpdf(30,:),'LineWidth',3);y1=get(gca,'ylim');hold on;plot([data_extended(31,1) data_extended(31,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=30','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,7)
plot(-60:60,kpdf(35,:),'LineWidth',3);y1=get(gca,'ylim');hold on;plot([data_extended(36,1) data_extended(36,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=35','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,8)
plot(-60:60,kpdf(40,:),'LineWidth',3);y1=get(gca,'ylim');hold on;plot([data_extended(41,1) data_extended(41,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=40','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);
subplot(3,3,9)
plot(-60:60,kpdf(45,:),'LineWidth',3);y1=get(gca,'ylim');hold on;plot([data_extended(46,1) data_extended(46,1)],[y1(1) y1(2)+0.1],'--', 'LineWidth',3,'MarkerSize',5);set(gca,'FontSize',20);
title('K=45','FontSize', 24)
xlabel('X','FontSize', 24);ylabel('PDF','FontSize', 24);

%%  EKF Evolution of Posterior PDF
figure
patchColor = [0 0.4470 0.7410];
[X,Y] = meshgrid(1:2:50,-60:60);
hFill = fill3(X,Y, kpdf(1:2:end,:), patchColor, 'LineWidth', 1, 'EdgeColor', patchColor,'FaceAlpha', 0.5);
grid on
hold on
plot(1:2:50,data_extended(2:2:end,1),'r','LineWidth',4)
title('EKF Evolution of Posterior PDF')
xlabel('K','FontSize', 24);ylabel('X','FontSize', 24);zlabel('PDF','FontSize', 24);set(gca,'FontSize',20);ylim([-60 60])

%% EKF True & Estimated State

figure
plot(1:50,data_extended(2:end,1),1:50,data_extended(2:end,2),'--or','LineWidth',3,'MarkerSize',10);set(gca,'FontSize',20);
title('EKF TRUE AND ESTIMATED STATE')
xlabel('K','FontSize', 24);ylabel('X','FontSize', 24);
legend('Truth','Estimate')

%% EKF Estimation Errorand 2-Sigma Bounds

figure 
sgtitle('EKF ESTIMATION ERROR AND 2-SIGMA BOUNDS')
plot(1:50,data_extended(2:end,1)-data_extended(2:end,2),'LineWidth',1.5);set(gca,'FontSize',20);hold on;
plot(1:50, 2.*sqrt(data_extended(2:end,3)),'r',1:50,-2.*sqrt(data_extended(2:end,3)),'r','LineWidth',1.5)
xlabel('K','FontSize', 24);ylabel('ERROR','FontSize', 24);
legend('Error','Bound')

%% Comparison
figure
sgtitle('Particle Filter - TRUE AND ESTIMATED STATE','FontSize', 24)
plot(1:50,data_extended(2:end,1),'b', 'LineWidth',3); hold on 
plot(1:50,data_particle(2:end,1),'--or','LineWidth',3)
legend('Truth','Estimate')
xlabel('K','FontSize', 24);ylabel('X','FontSize', 24);

%% Plot covariance of EKF
ti(3:2:51*2)=1:50;ti(2:2:51*2-1)=1:50;
P(1:2:101)=data_extended(:,3);P(2:2:101)=data_extended(2:end,4);
figure
plot(ti,P,0:50,data_extended(:,3),'*',1:50,data_extended(2:end,4),'o','LineWidth',1.5,'MarkerSize',10)
title('Covariance Matrix','FontSize', 24)
xlabel('K','FontSize', 24);ylabel('P','FontSize', 24);
legend('P ','P^+','P^-','FontSize', 24)
axes('Position',[.6 .4 .2 .2])
box on
plot(ti(50:60),P(50:60),26:30,data_extended(26:30,3),'*',25:30,data_extended(26:31,4),'o','LineWidth',1.5,'MarkerSize',5)