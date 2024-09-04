%% Define parameters
clear all, close all
global J1 J2 l1 l2 b1 b2 b1s m2 g Kinf;
J1 = 0.000027 + 0.00257; J2 = 4.15e-5;
l1 = 0.16; l2 = 0.075;
m2 = 0.0206; g = 9.81;
b1 = 0.01;  b2 = 0.001; b1s = 0.001;  % FRICTION PARAM: 0.01 0.001 0.02 

%% Define the linearized state-space model
gam = J1*J2 + J1*m2*l2^2 + J2*m2*l1^2;
a32 = m2^2*l1*l2^2*g / gam; a42 = m2*l2*g*(J1 + m2*l1^2) / gam;
a33 = -b1*(J2+m2*l2^2)/gam; a43 = -b1*m2*l1*l2/gam;
a34 = -b2*m2*l1*l2/gam;      a44 = -b2*(J1+m2*l1^2)/gam;
b31 = (J2 + m2*l2^2) / gam; b32 = m2*l1*l2 / gam;
b41 = m2*l1*l2 / gam; b42 = (J1 + m2*l1^2) / gam;

% Linearized state space model
A = [ 0 0 1 0;
      0 0 0 1;
      0 a32 a33 a34;
      0 a42 a43 a44 ];
B = [ 0 0 ;
      0 0 ;
      b31 b32;
      b41 b42 ];
x0 = [0; 5*pi/180; 0; 0];  % initial condition 

%% simulation parameters
x_eq = [0; 0; 0; 0];  % linearization equilibrium point
tf = 10; % simulation final time

%% Design infinte horizon LQR problem and solve
Q = diag([10, 1, 0, 0]);
R=[10 50 200];

tf1 = round(tf);
z = zeros(1, tf1+1);
tz=linspace(0, tf1, tf1+1);
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

for i=1:length(R)
    % compute the ARE to obtain retration's matrix Kinf
    [Kinf, Pinfinito, lam] = lqr(A, B(:,1), Q, R(i));

    % Simulate the feedback control with the nonlinear plant
    opts = odeset('MaxStep', 0.001, 'RelTol', 1e-8);
    [t,x_rad] = ode45(@sysCLode, [0 tf], x0, opts); 
    x_rad(:,2) = x_rad(:,2) + x_eq(2);
    x = x_rad .* 180./pi;

    %plot
    figure(1);
    subplot(2,2,1);
    plot(t,x(:,1));
    title('\theta_1');
    labels=arrayfun(@(x) sprintf('r = %d', x), R, 'UniformOutput', false);
    legend(labels);
    hold on;
    subplot(2,2,2);
    plot(t,x(:,2));
    title('\theta_2');
    legend(labels);
    hold on;
    subplot(2,2,3);
    plot(t,x(:,3));
    title('\theta_1^''');
    legend(labels);
    hold on;
    subplot(2,2,4);
    plot(t,x(:,4));
    title('\theta_2^''');
    legend(labels);
    hold on;
    
    u=-x*reshape(Kinf,4,1);
    figure(2);
    plot(t,u);
    title('Control law: u=\tau_1');
    legend(labels);
    hold on;
end





