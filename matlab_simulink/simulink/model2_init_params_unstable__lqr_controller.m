%% unstable model parameters
% M matrix -> it's symmetric
m1 = 96.8;
m2 = -3.57;
m3 = -3.57;
m4 = 0.258;

% C matrix
c1 = 0;
c2 = -50.8;
c3 = 0.436;
c4 = 2.20;

% K0 matrix
k11 = -901.0;
k12 = 35.17;
k13 = 35.17;
k14 = -12.03;

% K2 matrix
k21 = 0;
k22 = -87.06;
k23 = 0;
k24 = 3.50;

V = 5;
a1 = k11+k21*V^2;
a2 = k12+k22*V^2;
a3 = k13+k23*V^2;
a4 = k14+k24*V^2;

A1 = m2/(m2*m3-m1*m4);
A2 = m4/m2;
A3 = m1/m2;


%% state space
A = [
    0                                     0                             1                               0;
    0                                     0                             0                               1;
    A1*A2*a1-A1*a3                     A1*A2*a2-A1*a4                 A1*A2*c1*V-A1*c3*V                  A1*A2*c2*V-A1*c4*V;
    -A1*A2*A3*a1+A1*A3*a3-a1/m2     -A1*A2*A3*a2+A1*A3*a4-a2/m2     -A1*A2*A3*c1*V+A1*A3*c3*V-c1*V/m2     -A1*A2*A3*c2*V+A1*A3*c4*V-c2*V/m2;
];

B = [
    0;
    0;
    A1;
    -A3*A1
];

C = [1 0 0 0];
D = 0;

%% CCF form
[tf_num, tf_denum] = ss2tf(A,B,C,D);
A_ccf = [
    0                  1               0               0;
    0                  0               1               0;
    0                  0               0               1;
    -tf_denum(5)    -tf_denum(4)    -tf_denum(3)    -tf_denum(2)
];
B_ccf = [
    0;
    0;
    0;
    1
];

C_ccf = [tf_num(5)  tf_num(4)   tf_num(3)   tf_num(2)];

D_ccf = 0;


%% K matrix using LQR
Q = [
    1/0.01  0   0   0;
    0   1/0.04  0   0;
    0   0   0   0;
    0   0   0   0;
];
R = 1/4;

[K, S,P] = lqr(A_ccf, B_ccf, Q, R);

%% observer on the CCF form
observer_poles = [-26, -26, -26, -26]; % two times faster than the fastest system pole after using LQR State FB
L_t = acker(A_ccf', C_ccf', observer_poles);
L = L_t';

A_obs = A_ccf-L*C_ccf;
B_obs = [B_ccf L]; % input is [u Y]'
C_obs = eye(length(A_obs)); % since we want all the states
D_obs = zeros(size(B_obs));
