clc; clear;
% BlueROV2 Heavy Configuration Specifications 

g = 9.82;    % m/s^2
rho = 1000;  % kg/m3

m = 13.5;    % kg
V = 0.0134;  % m^3

Ix = 0.26;  % kg*m2
Iy = 0.23;  % kg*m2
Iz = 0.37;  % kg*m2

g_center =  [0, 0, 0];      % m
b_center =  [0, 0, -0.01];  % m

Xu_l = 13.7;   % Ns/m (Linear drag surge)
Yv_l = 0;      % Ns/m (Linear drag sway)
Zw_l = 33.0;   % Ns/m (Linear drag heave)
Kp_l = 0;      % Ns/m (Linear drag roll)
Mq_l = 0.8;    % Ns/m (Linear drag pitch)
Nr_l = 0;      % Ns/m (Linear drag yaw)

Xu_q = 141.0;  % Ns^2/m^2 (Quadratic drag surge)
Yv_q = 217.0;  % Ns^2/m^2 (Quadratic drag sway)
Zw_q = 190.0;  % Ns^2/m^2 (Quadratic drag heave)
Kp_q = 1.19;   % Ns^2/m^2 (Quadratic drag roll)
Mq_q = 0.47;   % Ns^2/m^2 (Quadratic drag pitch)
Nr_q = 1.5;    % Ns^2/m^2 (Quadratic drag yaw)

Xu_a = 6.36;   % kg (Added mass surge)
Yv_a = 7.12;   % kg (Added mass sway)
Zw_a = 18.68;  % kg (Added mass heave)
Kp_a = 0.189;  % kg*m^2 (Added mass roll)
Mq_a = 0.135;  % kg*m^2 (Added mass pitch)
Nr_a = 0.222;  % kg*m^2 (Added mass yaw)


% ------------------------------Dynamics-----------------------------------


M_rb = eye(6).*[m, m, m, Ix, Iy, Iz];                % Rigid body mass
M_a = eye(6).*[Xu_a, Yv_a, Zw_a, Kp_a, Mq_a, Nr_a];  % Added mass
M = M_rb + M_a;  % Mass Matrix

D_l = eye(6).*[Xu_l, Yv_l, Zw_l, Kp_l, Mq_l, Nr_l];  % Linear drag matrix
D_q = eye(6).*[Xu_q, Yv_q, Zw_q, Kp_q, Mq_q, Nr_q];  % Quadratic drag matrix
D = D_l + D_q;  % Hydrodynamic damping

W = m*g;      % Weight
B = rho*g*V;  % Buoyancy

%    x  y  z  f  t  r
p = [0, 0, 0, 0, 0, 0];  % SİL BENİ

x_g = g_center(1);
y_g = g_center(2);
z_g = g_center(3);
x_b = b_center(1);
y_b = b_center(2);
z_b = b_center(3);

% Restoring force
g_n = [                    (W-B)*sin(p(5));
                     -(W-B)*cos(p(5))*sin(p(4));
                     -(W-B)*cos(p(5))*cos(p(4));
  -(y_g*W-y_b*B)*cos(p(5))*cos(p(4)) + (z_g*W-z_b*B)*cos(p(5))*sin(p(4));
        (z_g*W-z_b*B)*sin(p(5)) + (x_g*W-x_b*B)*cos(p(5))*cos(p(4));
       -(x_g*W-x_b*B)*cos(p(5))*sin(p(4)) - (y_g*W-y_b*B)*sin(p(5)) ];

u = [0, 0, 1, 0, 0, 0]';  % SİL BENİ

% Thruster allocation matrix
K = [ 0.7071,  0.7071, -0.7071, -0.7071,       0,       0,       0,       0;
     -0.7071,  0.7071, -0.7071,  0.7071,       0,       0,       0,       0;
           0,       0,       0,       0, -1.0000, -1.0000, -1.0000, -1.0000;
           0,       0,       0,       0,  0.2180,  0.2180, -0.2180, -0.2180;
           0,       0,       0,       0,  0.1200, -0.1200,  0.1200, -0.1200;
     -0.1888,  0.1888,  0.1888, -0.1888,       0,       0,       0,       0 ];

% Control allocation matrix
A = [-1,  1,  0,  0,  0,  1;
     -1, -1,  0,  0,  0, -1;
      1,  1,  0,  0,  0, -1;
      1, -1,  0,  0,  0,  1;
      0,  0, -1, -1,  1,  0;
      0,  0, -1,  1,  1,  0;
      0,  0, -1,  1, -1,  0;
      0,  0, -1, -1, -1,  0 ];

a = A*u;                              % Allocated control signal

t = (80 ./ ( 1 + exp(-4.*a.^3) )) - 40;  % Thrust
K
t
To = K*t;  % Force and moment from propellers
To

Td = [0, 0, 0, 0, 0, 0]';  % SİL BENİ
v = [0, 0, 0, 0, 0, 0]';  % SİL BENİ

v_dot = M^(-1) * (To + Td - D*v - g_n);
v_dot

% -----------------------------Kinematics----------------------------------


% Transformation matrix of linear velocity
R_n = zeros([3, 3]);

R_n(1, 1) = cos(p(4))*cos(p(5));
R_n(1, 2) = -sin(p(6))*cos(p(4)) + cos(p(6))*sin(p(5))*sin(p(4));
R_n(1, 3) = sin(p(6))*sin(p(4)) + cos(p(6))*cos(p(4))*sin(p(5));

R_n(2, 1) = sin(p(4))*cos(p(5));
R_n(2, 2) = cos(p(6))*cos(p(4)) + sin(p(4))*sin(p(5))*sin(p(6));
R_n(2, 3) = -cos(p(6))*sin(p(4)) + sin(p(5))*sin(p(6))*cos(p(4));

R_n(3, 1) = -sin(p(5));
R_n(3, 2) = cos(p(5))*sin(p(4));
R_n(3, 3) = cos(p(5))*cos(p(4));

% Transformation matrix of angular velocity
T_n = zeros([3, 3]);

T_n(1, 1) = 1;
T_n(1, 2) = sin(p(6))*sin(p(5))/cos(p(5));
T_n(1, 3) = cos(p(4))*sin(p(5))/cos(p(5));

T_n(2, 1) = 0;
T_n(2, 2) = cos(p(4));
T_n(2, 3) = sin(p(4));

T_n(3, 1) = 0;
T_n(3, 2) = sin(p(4))/cos(p(5));
T_n(3, 3) = cos(p(4))/cos(p(5));

J_n = [R_n, zeros([3,3]); zeros([3,3]), T_n];  % Transformation to the body velocity

p_dot = J_n * v;