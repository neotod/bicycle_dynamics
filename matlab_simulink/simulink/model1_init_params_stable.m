%% stable model parameters
b = 1.02;
c = 0.08;
h = 0.8603;
lambda = pi/10;
V = 5.7; % m/s
a = 0.347;
g = 9.81;
m = 84;

k1 = b^2/((V^2*sin(lambda) - b*g*cos(lambda))*m*a*c*sin(lambda));
k2 = b*g/(V^2*sin(lambda) - b*g*cos(lambda));

%% state space -> ccf
A = [
    0                                                       1;
    g/h - k2*(V^2*h-a*c*g)*sin(lambda)/(b*h^2)    -a*V*sin(lambda)/(b*h)*k2
];
B = [
    0;
    1
];

C = [(V^2*h-a*c*g)*sin(lambda)/(b*h^2)*k1 a*V*sin(lambda)/(b*h)*k1];
D = 0;

%% observer

observer_poles = [-12, -12];
L_t = acker(A', C', observer_poles);
L = L_t';
A_obs = A-L*C;
B_obs = [B L]; % input is [u Y]'
C_obs = eye(length(A_obs)); % since we want all the states
D_obs = zeros(size(B_obs));


%A_obs = [
  %-15.0970    0.0719
 %-133.4168   -8.9030
%];

%B_obs = [
         %0    0.4831
    %1.0000   -2.2958
%];

%C_obs = [
     %1     0;
     %0     1
%];
%D_obs = [
     %0     0
     %0     0
%];