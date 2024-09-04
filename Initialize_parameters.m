%main for Furuta's pendulum by simulink
clear all, close all

%% Define parameters
global J1 J2 l1 l2 m2 b1 b2 b1s g Kinf;
%mechanical parameters
J1 = 0.000027+0.00257; J2 = 4.15e-5;
l1 = 0.16; l2 = 0.075;
m2 = 0.0206; g = 9.81;
b1 = 0.01;  b2 = 0.001; b1s = 0.00;  % FRICTION PARAM: 0.01 0.001 0.02

rad2deg = 180/pi;   deg2rad = pi/180;

Ki = 2; %current-loop voltage gain [A/V]
Kt = 0.071; %torque constant [Nm/A]
Vmax = 3; %current-loop saturation of the reference [V]
b_DAC = 14; %number of bit of the digital controller
OutRange = 20; %output range of the controller [V]
q_DAC = OutRange/2^b_DAC; %quantization step of the controller
encoder_res = 2*pi/2000; %resolution of the encoders

x0 = [0; 5*deg2rad; 0; 0]; %pendulum initial condition

Td=0.01; %speed filter parameters
epsilon=1/sqrt(2); %speed filter parameters
ts_ZH=0.001; %zero order hold sample time

%pulse disturbance parameters
tr = 3; %torque disturbance t start [s]
tau2 = 0.00; %torque disturbance magnitude [Nm]
%% Define the linearized state-space model
gam = J1*J2 + J1*m2*l2^2 + J2*m2*l1^2;
a32 = m2^2*l1*l2^2*g / gam; a42 = m2*l2*g*(J1 + m2*l1^2) / gam;

a33 = -b1*(J2+m2*l2^2)/gam; a43 = -b1*m2*l1*l2/gam;
a34 = -b2*m2*l1*l2/gam;      a44 = -b2*(J1+m2*l1^2)/gam;

b31 = (J2 + m2*l2^2) / gam; b32 = m2*l1*l2 / gam;
b41 = m2*l1*l2 / gam; b42 = (J1 + m2*l1^2) / gam;
%{
A = [ 0 0 1 0;
      0 0 0 1;
      0 a32 0 0;
      0 a42 0 0]; 
%}
A = [ 0 0 1 0;
      0 0 0 1;
      0 a32 a33 a34;
      0 a42 a43 a44];
B = [ 0 0 ;
      0 0 ;
      b31 b32;
      b41 b42 ];

%Rescaling matrix for B, to derive voltage control
T = [Kt*Ki , 0 ;
       0   , 1 ];
BT = B*T;
ident = eye(4);
C = ident;
invC = inv(C);
D = zeros(4,2);
x_eq = [0; 0; 0; 0];

%% Design the infinte horizon LQR
% Bryson rule
alfaM = 30*deg2rad;
betaM = 10*deg2rad;
a_dotM = alfaM*2*pi*10;
b_dotM = betaM*2*pi*10;
U_Max = 0.5*Ki*Kt;

%Q = diag([1/alfaM^2, 1/betaM^2, 1/a_dotM^2, 1/b_dotM^2]);
%r = 1/U_Max^2;
Q = diag([3.6476   32.8281    0.0009    0.0083]); %[1 1 100 1]
r = 198.3733; %10

[Kinf, Pinf, lam] = lqr(A, BT(:,1), Q, r);

%% pole placement alternative
% alternatively, one could directly allocate the poles in arbitrary positions
% Here we allocate poles while varying w in different simulations

s=tf('s');
w=7; %0.1 1 10 100 ?
allocated_poles = [-w*sqrt(2)/2+1i*w*sqrt(2)/2 ,
                   -w*sqrt(2)/2-1i*w*sqrt(2)/2 ,
                   -w*sqrt(3)/2+1i*w*1/2       ,
                   -w*sqrt(3)/2-1i*w*1/2];
K = place( A, BT(:,1), allocated_poles ); %feedback matrix
%to be used to compare the behaviour with Kinf(LQR), while varying w

eigenval = eig( A - BT(:,1)*K );
Co=ctrb( A, B(:,1) );
unco = length(A) - rank(Co);

