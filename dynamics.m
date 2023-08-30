%% dynamics

% function [t_dot, w_dot] = dYNAMICS(t, t1,t2,t3,t4, w1,w2,w3)

function output = dynamics(t,states)

% q = [q1 q2 q3 q4]';     % Current Attitude quaternions
% w = [w1 w2 w3]';

t1 = states(1); t2 = states(2);
t3 = states(3); t4 = states(4);
w1 = states(5); w2 = states(6);
w3 = states(7);

% Lyapunov based Controller
J = [10 -3 -7; -3 18 2; -7 2 8];
wn = 0.1;
zeta = 0.7;
Kp = (wn^2)*J;
Kd = 2*zeta*wn*J;

tgq = [t1 t2 t3 t4]';       % To-go quaternions

%*********************************** 
% Output = d1,d2,d3,d4 & s1,s2,s3,s4
% Desired attitude quaternions
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

s1_dot = (d1_ddot*d4 + d1_dot*d4_dot) + (d2_ddot*d3 + d2_dot*d3_dot) ...
            - (d3_ddot*d2 + d3_dot*d2_dot) - (d4_ddot*d1 + d4_dot*d1_dot);
s2_dot = -(d1_ddot*d3 + d1_dot*d3_dot) + (d2_ddot*d4 + d2_dot*d4_dot) ...
            + (d3_ddot*d1 + d3_dot*d1_dot) - (d4_ddot*d2 + d4_dot*d2_dot);
s3_dot = (d1_ddot*d2 + d1_dot*d2_dot) - (d2_ddot*d1 + d2_dot*d1_dot) ...
            + (d3_ddot*d4 + d3_dot*d4_dot) - (d4_ddot*d3 + d4_dot*d3_dot);
s4_dot = 0;

% S = [s1 s2 s3 s4]';
% S_dot = [s1_dot s2_dot s3_dot s4_dot]';

%********************************

sx = [0     s3      -s2     s1;
     -s3    0       s1      s2;
     s2     -s1     0       s3;
     -s1    -s2     -s3     0];

omega = [0 w3 -w2 -w1;
        -w3 0 w1 -w2;
        w2 -w1 0 -w3;
        w1 w2 w3 0];

wx = [0 -w3 w2; w3 0 -w1; -w2 w1 0];
w = [w1 w2 w3]';
T = [t1, t2, t3]';
S = [s1 s2 s3]';
S_dot = [s1_dot s2_dot s3_dot]';

% % CONTROL LAW (NON-LINEAR)
u1 = wx*J*w + Kp*T - Kd*w + 2*( Kd*S + J*S_dot);

% % CONTROL LAW (Classical PD controller)
% u2 =  wx*J*w + Kp*T - Kd*w;

%% % FINAL DYNAMICS

% t_dot = [t1_dot; t2_dot; t3_dot; t4_dot]
% w_dot = inv(J)*(-wx*J*w + u)

w_dot = J\(-wx*J*w + u1);           % Non-linear controller
% w_dot = J\(-wx*J*w + u2);         % Classical PD controller
t_dot = ( sx + s4*eye(4) + 0.5*omega ) * tgq;

t_dot = reshape(t_dot, [], 1);
w_dot = reshape(w_dot, [], 1);

output = [t_dot; w_dot]

D = [-d4 -d3 d2 d1; d3 -d4 -d1 d2;  -d2 d1 -d4 d3;  d1 d2 d3 d4];  
new_q = D\[t1,t2,t3,t4]';

e_q1 = new_q(1)-d1;
e_q2 = new_q(2)-d2;
e_q3 = new_q(3)-d3;
e_q4 = new_q(4)-d4;

% % q plot 1,2
plot(t,new_q(1),'.-r', t,new_q(2),'.-b', t,new_q(3),'.-g', t,new_q(4),'.-k')
hold on

% % Control plot 1
% plot(t,u1(1),'.-r', t,u1(2),'.-b', t,u1(3),'.-k')
% hold on

% % Error plot 1,2
% plot(t,e_q1,'.-r', t,e_q2,'.-b', t,e_q3,'.-g', t,e_q4,'.-k')
% hold on

% % Control plot 2
% plot(t,u2(1),'.-r', t,u2(2),'.-b', t,u2(3),'.-k')
% hold on

end