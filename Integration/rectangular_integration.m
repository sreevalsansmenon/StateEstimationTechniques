function integral=rectangular_integration(funs,time,stepsize,x_t)
% Program by Sreevalsan S Menon(sm2hm@mst.edu)
% Rectangular Integration
% This program can be used to perform rectangular integration the 'funs' is
% the function that has to be integrated, it should be entered as a
% function handle. eg @(t) cos(t), @(t)[cos(t);sin(t)] @(t,x1,x2)[x1;x2]
% Note t should always be provided in function handle as first variable.
% 'time' is the total time of integration.'stepsize' is the integration
% step size.'x_0' is the initial conditions.
n=1;                                                                    % Variable to point step number integral values 
integral(:,n)=x_t;                                                      % Variable to save integral value
for ti=0:stepsize:time-stepsize                                         % Loop to find integral
    n=n+1;                                                              % Increment step number
    brace=num2cell([ti ;integral(:,n-1)]);                              % Converting numbers to cell
    integral(:,n)=integral(:,n-1)+funs(brace{1:nargin(funs)})*stepsize; % Finding integral for ti
end                                                                     % Loop Terminates
