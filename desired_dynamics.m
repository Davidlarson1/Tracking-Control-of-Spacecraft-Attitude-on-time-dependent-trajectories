%%
% t = 0:0.0005:100;
% close all
% clear all
% clc
t_collect = [];
d1_collect = [];
d2_collect = [];
d3_collect = [];
d4_collect = [];
% for i = 1:0.0005:100.0005

for i = 1:0.5:100.5
    t = i-1;

lmbd1 = 1/sqrt(14);
lmbd2 = 2/sqrt(14);
lmbd3 = 3/sqrt(14);

a = 0;
b = 0;
c = 2*pi*(0.0001);
e = -(4/3)*pi*(0.000001);

alpha = a + b*t + c*t^2 + e*t^3;
alpha_dot = b + 2*c*t + 3*e*t^2;
alpha_ddot = 2*c + 6*e*t;

d1 = lmbd1*sin(alpha/2);
d2 = lmbd2*sin(alpha/2);
d3 = lmbd3*sin(alpha/2);
d4 = cos(alpha/2);

% D = [-d4 -d3 d2 d1; d3 -d4 -d1 d2;  -d2 d1 -d4 d3;  d1 d2 d3 d4];  

d1_dot = (lmbd1/2)*alpha_dot*cos(alpha/2);
d2_dot = (lmbd2/2)*alpha_dot*cos(alpha/2);
d3_dot = (lmbd3/2)*alpha_dot*cos(alpha/2);
d4_dot = -(alpha_dot/2)*sin(alpha/2);

d1_ddot = (lmbd1/2)*( (alpha_ddot*cos(alpha/2)) - ((alpha_dot^2)/2)*sin(alpha/2) );
d2_ddot = (lmbd2/2)*( (alpha_ddot*cos(alpha/2)) - ((alpha_dot^2)/2)*sin(alpha/2) );
d3_ddot = (lmbd3/2)*( (alpha_ddot*cos(alpha/2)) - ((alpha_dot^2)/2)*sin(alpha/2) );
d4_ddot = -(0.5*alpha_ddot*sin(alpha/2)) - (0.25*(alpha_dot^2)*cos(alpha/2));

s1 = d1_dot*d4 + d2_dot*d3 - d3_dot*d2 - d4_dot*d1;
s2 = -d1_dot*d3 + d2_dot*d4 + d3_dot*d1 - d4_dot*d2;
s3 = d1_dot*d2 - d2_dot*d1 + d3_dot*d4 - d4_dot*d3;
s4 = d1_dot*d1 + d2_dot*d2 + d3_dot*d3 + d4_dot*d4;
s4 = 0;

t_collect = [t_collect,t]
d1_collect = [d1_collect, d1];
d2_collect = [d2_collect, d2];
d3_collect = [d3_collect, d3];
d4_collect = [d4_collect, d4];
end

figure(1)
plot(t_collect,d1_collect,'r')
hold on
plot(t_collect,d2_collect,'g')
hold on
plot(t_collect,d3_collect,'b')
hold on
plot(t_collect,d4_collect,'k')
hold off
title('Desired Attitude as a function of time')
xlabel('Time [sec]')
ylabel('Desired Attitude')
legend('d1','d2','d3','d4')