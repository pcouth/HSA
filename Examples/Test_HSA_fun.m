%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
%
% A script to test some output functions from HSA generation using the 
% Test_HSA.m file 
%--------------------------------------------------------------------------
%% 1D1P
clc; close all

p = -0.5;
tf = 10;
x0 = .7;
[t,x] = ode45(@(t,x) SN_gen(x,[],p),[0 tf],x0);

figure(1)
plot(t,x,'b')

%% 1DMP
clc; close all

p = [0 -0.2 1 1];
tf = 10;
x0 = -3;
[t,x] = ode45(@(t,x) ODMP_fp(x,[],p),[0 tf],x0);

figure(1)
plot(t,x,'b')

%% MD1P
clc; close all

p = -0.5;
tf = 10;
x0 = [-3 0];
[t,x] = ode45(@(t,x) MD1P_fp(x,[],p),[0 tf],x0);

figure(1)
plot3(t,x(:,1),x(:,2),'b')
xlabel('Time'); ylabel('x_1'); zlabel('x_2')
grid on

%% MDMP
clc; close all

p = [1 -0.5 0 1];
% p = [-0.5 0 1 0];
tf = 10;
x0 = [-1 0];
% [t,x] = ode45(@(t,x) MDMP_gen(x,[],p),[0 tf],x0);
[t,x] = ode45(@(t,x) MDMP_fp(x,[],p),[0 tf],x0);

figure(1)
plot3(t,x(:,1),x(:,2),'b')
xlabel('Time'); ylabel('x_1'); zlabel('x_2')
grid on


