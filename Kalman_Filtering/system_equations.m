% Program by Sreevalsan S Menon(sm2hm@mst.edu)

function dxdt =system_equations(t,x0)
% System Dynamics
rho_0=0.0034;
g=32.2;
k=22000;

x1_dot=x0(2);
x2_dot=rho_0*exp(-x0(1)/k)*x0(2)^2/(2*x0(3)) - g;
x3_dot=0; 
dxdt=[x1_dot; x2_dot; x3_dot];
