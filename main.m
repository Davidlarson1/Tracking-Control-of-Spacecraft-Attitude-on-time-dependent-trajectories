%% MAIN
close all
clear all
clc

w_initial = [0.1 0.2 0.1]';        % Initial conditions

% ****************************************
% Desired attitude quaternions
t = 0;
lmbd1 = 1/sqrt(14);
lmbd2 = 2/sqrt(14);
lmbd3 = 3/sqrt(14);

a = 0;
b = 0;
c = 2*pi*(0.0001);
e = -(4/3)*pi*(0.000001);

alpha = a + b*t + c*t^2 + e*t^3;
% alpha_dot = b + 2*c*t + 3*e*t^2;
% alpha_ddot = 2*c + 6*e*t;

d1 = lmbd1*sin(alpha/2);
d2 = lmbd2*sin(alpha/2);
d3 = lmbd3*sin(alpha/2);
d4 = cos(alpha/2);
%**********************************

q_initial = [0 0 0 1]';
% d_initial = [d1 d2 d3 d4];

D = [-d4 -d3 d2 d1; d3 -d4 -d1 d2;  -d2 d1 -d4 d3;  d1 d2 d3 d4];  
t_initial = D*q_initial;

w_0 = w_initial;
w1_0 = w_0(1);  w2_0 = w_0(2);  w3_0 = w_0(3);
tg_0 = t_initial;
t1_0 = tg_0(1); t2_0 = tg_0(2); t3_0 = tg_0(3); t4_0 = tg_0(4);

states_0 = [t1_0,t2_0,t3_0,t4_0, w1_0,w2_0,w3_0];

% ODE45 INTEGRATION
[t,y] = ode45(@dynamics, [0 100], states_0)

t1 = y(:,1);
t2 = y(:,2);
t3 = y(:,3);
t4 = y(:,4);

w1 = y(:,5);
w2 = y(:,6);
w3 = y(:,7);

figure()
plot(t,t1,'r')
hold on
plot(t,t2,'b')
hold on
plot(t,t3,'g')
hold on
plot(t,t4,'k')
hold off
title('To-go quaternions variation with time')
xlabel('Time [s]')
ylabel('To-go quaternions')
legend('t_1','t_2','t_3','t_4')

figure()
plot(t,w1,'r')
hold on
plot(t,w2,'b')
hold on
plot(t,w3,'g')
hold off
title('Angular rates of Satellite')
xlabel('Time [s]')
ylabel('Angular rates (\omega)')
legend('\omega_1','\omega_2','\omega_3')
