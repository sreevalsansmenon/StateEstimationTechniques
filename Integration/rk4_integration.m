function integral=rk4_integration(funs,time,stepsize,x_0)
% Program by Sreevalsan S Menon(sm2hm@mst.edu)
% Runga Kutta 4 Integration
% This program can be used to perform rk4 integration the 'funs' is
% the function that has to be integrated, it should be entered as a
% function handle. eg @(t) cos(t), @(t)[cos(t);sin(t)] @(t,x1,x2)[x1;x2]
% Note t should always be provided in function handle as first variable.
% 'time' is the total time of integration.'stepsize' is the integration
% step size.'x_0' is the initial conditions.
n=1;                                                        % Variable to point step number integral values 
integral(:,n)=x_0;                                          % Variable to save integral value
for ti=0:stepsize:time-stepsize                             % Loop to find integral
    n=n+1;                                                  % Increment step number
    brace=num2cell([ti;integral(:,n-1)]);                   % Converting numbers to cell
    dx_1=funs(brace{1:nargin(funs)})*stepsize;              % Finding delta x1
    brace=num2cell([ti+0.5*stepsize ;integral(:,n-1)+0.5*dx_1]);     % Converting numbers to cell
    dx_2=funs(brace{1:nargin(funs)})*stepsize;              % Finding delta x2
    brace=num2cell([ti+0.5*stepsize ;integral(:,n-1)+0.5*dx_2]);     % Converting numbers to cell
    dx_3=funs(brace{1:nargin(funs)})*stepsize;              % Finding delta x3
    brace=num2cell([ti+stepsize; integral(:,n-1)+dx_3]);         % Converting numbers to cell
    dx_4=funs(brace{1:nargin(funs)})*stepsize;              % Finding delta x4
    integral(:,n)=integral(:,n-1)+(dx_1+2*dx_2+2*dx_3+dx_4)/6;%  Finding integral for ti
end                                                         % Loop Terminates