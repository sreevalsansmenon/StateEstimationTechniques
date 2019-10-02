% Program by Sreevalsan S Menon(sm2hm@mst.edu)

clc                     % Clear command window
clear                   % clear workspace

%  Use different step size size and integration methods for cos(t) and find error percentage for time 0-1 seconds

stepsize=[0.1;0.05;0.025];                  % Step size for comparison
for i=1:size(stepsize,1)                    % Loop in step size
    rectangularintegral=rectangular_integration(@(t)cos(t),1,stepsize(i),0);    % Finding integral using rectangular integration
    trapezoidalintegral=trapezoidal_integration(@(t)cos(t),1,stepsize(i),0);    % Finding integral using trapezoidal integration
    rk4integral=rk4_integration(@(t)cos(t),1,stepsize(i),0);                    % Finding integral using rk4 method
    fprintf('Percentage error in rectangular integration for %d stepsize is %d\n',stepsize(i),100*abs((sin(1)-rectangularintegral(end))/sin(1) )) % Print percentage error 
    fprintf('Percentage error in trapezoidal integration for %d stepsize is %d\n',stepsize(i),100*abs((sin(1)-trapezoidalintegral(end))/sin(1) )) % Print percentage error
    fprintf('Percentage error in runge-kutta integration for %d stepsize is %d\n\n',stepsize(i),100*abs((sin(1)-rk4integral(end))/sin(1) ))       % Print percentage error
end                                                                             % Loop terminates

%% Homework Part 2
% Integrating system of equations
% Finding closed form solution

syms x1(t) x2(t)                                                                % Define differntial variable
odes=[diff(x1) == x2;diff(x2) == -x1];                                          % Given system of equations
Initialconds = [x1(0) == 1;x2(0) == 2];                                         % initial conditions
[x1Sol(t), x2Sol(t)] = dsolve(odes,Initialconds);                               % Solving for closed form solutions
t=0:0.01:5;x1_t=eval(x1Sol(t));x2_t=eval(x2Sol(t));                             % Findin integral value using closed form solution

% Find integral using Numerical Integration
stepsize=[ 0.1 ;0.05 ;0.025];                           % Step size for comparison
for i=1:size(stepsize,1)                                % Loop in step size
    rectangularintegral=rectangular_integration(@(t,x1,x2)[x2;-x1],5,stepsize(i),[1;2]);    % Finding integral using rectangular integration
    trapezoidalintegral=trapezoidal_integration(@(t,x1,x2)[x2;-x1],5,stepsize(i),[1;2]);    % Finding integral using trapezoidal integration
    rk4integral=rk4_integration(@(t,x1,x2)[x2;-x1],5,stepsize(i),[1;2]);                    % Finding integral using rk4 method
    fprintf('Percentage error in x1(t) rectangular integration for %d stepsize is %d\n',stepsize(i),100*abs((x1_t(end)-rectangularintegral(1,end))/x1_t(end) )) % Print percentage error 
    fprintf('Percentage error in x1(t) trapezoidal integration for %d stepsize is %d\n',stepsize(i),100*abs((x1_t(end)-trapezoidalintegral(1,end))/x1_t(end) )) % Print percentage error
    fprintf('Percentage error in x1(t) runge-kutta integration for %d stepsize is %d\n',stepsize(i),100*abs((x1_t(end)-rk4integral(1,end))/x1_t(end) ))       % Print percentage error
    fprintf('Percentage error in x2(t) rectangular integration for %d stepsize is %d\n',stepsize(i),100*abs((x2_t(end)-rectangularintegral(2,end))/x2_t(end) )) % Print percentage error 
    fprintf('Percentage error in x2(t) trapezoidal integration for %d stepsize is %d\n',stepsize(i),100*abs((x2_t(end)-trapezoidalintegral(2,end))/x2_t(end) )) % Print percentage error
    fprintf('Percentage error in x2(t) runge-kutta integration for %d stepsize is %d\n\n',stepsize(i),100*abs((x2_t(end)-rk4integral(2,end))/x2_t(end) ))       % Print percentage error
    % Figure to plot x1(t)
    figure
    plot(t,x1_t)
    hold on
    plot(0:stepsize(i):5,rectangularintegral(1,:),'*');set(gca,'FontSize',30);
    plot(0:stepsize(i):5,trapezoidalintegral(1,:),'o')
    plot(0:stepsize(i):5,rk4integral(1,:),'x')
    xlabel('Time in seconds');ylabel('X1(t)')
    legend('Closed Form', 'Rectangular Integration','Trapezoidal Integration','Runge Kutta 4')
    title(['Numerical integration x1(t) using step size ' num2str(stepsize(i))])
    hold off
    
    % Figure to plot x2(t)
    figure
    plot(t,x2_t)
    hold on
    plot(0:stepsize(i):5,rectangularintegral(2,:),'*');set(gca,'FontSize',30);
    plot(0:stepsize(i):5,trapezoidalintegral(2,:),'o')
    plot(0:stepsize(i):5,rk4integral(2,:),'x')
    xlabel('Time in seconds');ylabel('X2(t)')
    legend('Closed Form', 'Rectangular Integration','Trapezoidal Integration','Runge Kutta 4')
    title(['Numerical integration x2(t) using step size ' num2str(stepsize(i))])
    hold off
end                                                                             % Loop terminates

