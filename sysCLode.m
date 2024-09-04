function dx = sysCLode(t, x)
global J1 J2 m2 l1 l2 b1 b2 b1s g Kinf;
 
th2 = x(2);
dth1 = x(3);
dth2 = x(4);

%% Define nonlinear model: M*ddq + C*dq + G = u + h
M = [J1+m2*l1^2 + m2*l2^2*sin(th2)^2, -m2*l1*l2*cos(th2);
      -m2*l1*l2*cos(th2)            ,     J2 + m2*l2^2 ];

C = [ m2*l2^2*dth2*sin(th2)*cos(th2)+b1, ...
        m2*l2^2*dth1*sin(th2)*cos(th2) + m2*l1*l2*dth2*sin(th2);
        -m2*l2^2*dth1*sin(th2)*cos(th2), +b2];
    
G = [0; -m2*l2*g*sin(th2)];

h = [-b1s*sign(dth1) ; 0];

%% Generate a pulse disturbance
if t>3 && t < 3.1
    tau2 = 0.01; %pulse disturbance
else
    tau2 = 0;
end

%% Define the feedback control & disturbance inputs
u = [-Kinf*x; tau2];

%% Update derivatives
dq = [dth1; dth2];
ddq = M\(-C*dq-G+u+h);
dx = [dq; ddq];
end