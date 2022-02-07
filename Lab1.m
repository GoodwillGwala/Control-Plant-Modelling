%%%%%%%%%%%%%%%%%%%%%%% Backward Euler Method %%%%%%%%%%%%%%%%%%%
%AUTHORS : Group 59
%LAB 1 & 2
%
%Acknowledgements: Our sincerest gratitude to Mubanga Mofu for his guidance
%during the course of the completion of this lab exercise.
clc;
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%% Variable Declaration %%%%%%%%%%

CL_file = 'CL.mat';
OL_file = 'OL.mat';
TD_file = 'TD.mat';
PD1_file = 'PD1.mat';
PD2_file = 'PD2.mat';

[CL] = load(CL_file); [OL] = load(OL_file); [TD] = load(TD_file); [PD1] = load(PD1_file); [PD2] = load(PD2_file);


[Variable] = load('Variable.mat');
%ensure that dt_sub is an integer else use dt = 0.2

Variable.x0 = [0;1];
Variable.X_OL( :,1) = Variable.x0;  %Open Loop System
Variable.X_GC( :,1) = Variable.x0; %define proportional gain

if ~isinf(Variable.m) && floor(Variable.m) == Variable.m
    dt = Variable.dt_sub;
else
    dt = Variable.tau;
    Variable.m = Variable.tau/dt;
end

%define the differential consant variables

P0 = 1/(dt); P1 = 1+8*dt; P2 = P0-((20*dt)/(P1));

A11 = P0/P2; A12 = (-30*dt)/(P1*P2); A13 = 1/(P1*P2);
B_1  = dt/(P1*P2);

A21 = (A11*P0)-P0; A22 = A12*P0; A23 = A13*P0;
B_2  = B_1*P0;

%Gain Scheduling

switch Variable.x0(1)
    case 0
        Kp = 1.1;
    otherwise
        Kp = 1;
end

for n=1:Variable.T/dt
    tF(n+1) = n*dt;
    if n-Variable.m+1<1 %ensure that the system is causal
        m_del = 1;
    else
        m_del = n-Variable.m+1;
    end
    %Closed Loop
    Variable.X(1,n+1) = A11*Variable.X(1,n) + A12*Variable.X(1,m_del) + A13*Variable.X(2,n) + B_1*Variable.u;
    Variable.X(2,n+1) = A21*Variable.X(1,n) + A22*Variable.X(1,m_del) + A23*Variable.X(2,n) + B_2*Variable.u;
    %Open Loop
    Variable.X_OL(n+1) = (1/(1-2*dt)) * Variable.X_OL(n) + (6*dt)/(1-2*dt)*Variable.u;
    %Gain Controller
    Variable.X_GC(1,n+1) = A11*Variable.X_GC(1,n) + Kp*A12*Variable.X_GC(1,m_del) + A13*Variable.X_GC(2,n) + B_1*Variable.u;
    Variable.X_GC(2,n+1) = A21*Variable.X_GC(1,n) + Kp*A22*Variable.X_GC(1,m_del) + A23*Variable.X_GC(2,n) + B_2*Variable.u;


end

%Closed Loop Output
Y = 60.*Variable.X(1,:) + 6.*Variable.X(2,:);
save ('CL_ML.mat', 'tF', 'Y');
%Open Loop Output
Y_OL = 2.*Variable.X_OL +6*Variable.u;
save ('OL_ML.mat', 'tF', 'Y_OL')
%Compensated KP Controller Output
Y_KP = Kp*(60.*Variable.X_GC(1,:) + 6.*Variable.X_GC(2,:));
save ('GC_ML.mat', 'tF', 'Y_KP')
Y_Y_KP = [Y ; Y_KP];

figure;
grid on;
y_thres = 200;
Steady_state = Y(end);
Room_Temperature = Y(1);
plot(tF,[Y;CL.Y]);
grid on;
line([0,tF(1,end)],[y_thres,y_thres],'Color','red','LineStyle','--');
line([0,tF(1,end)],[ Steady_state,Steady_state],'Color','green','LineStyle','--');
title('Essential Oil Extractor Temperature System Response')
xlabel('Time [s]')
ylabel('Temperature [°C]')
legend('System Response SIMULINK','System Response','Maximum Temperature', 'Steady State');
axis ([0 10 0 250])


% Backward Euler vs Gain Controller Plot

figure;

Steady_state = Y_KP(end);
plot(tF,[Y_Y_KP]);
grid on;
line([0,tF(1,end)],[200,200],'Color','red','LineStyle','--');
line([0,tF(1,end)],[ Steady_state,Steady_state],'Color','green','LineStyle','--');
line([0,tF(1,end)],[ 100,100],'Color','magenta','LineStyle','--');
title('Uncompensated vs Compensated System using Gain Controller x_1(3.34)')
xlabel('Time [s]')
ylabel('Temperature [°C]')
legend('System Response', 'Compensated System Response','Maximum Temperature', 'New Steady State', 'Minimum Temperature');
axis ([0 10 0 250])


figure;
grid on;
plot(tF,[Y_OL'; OL.Y]);
grid on;
line([0,tF(1,end)],[y_thres,y_thres],'Color','red','LineStyle','--');
title('Open Loop System Response')
xlabel('Time [s]')
ylabel('Temperature [°C]')
legend('System Response','System Response SIMULINK', 'Maximum Temperature');
axis ([0 100 0 1e14])


%PADE APPROX

p = 2/Variable.tau;

%1st Order Coeffs
a3 = 1; a2 = 8+p; a1 = 8*p-50; a0 = 10*p;

b3 = 0; b2 = 6; b1 = 60+6*p; b0 = 60*p;

%Transfer Function with Pade Approximation delay

b = [b3 b2 b1 b0];
a = [a3 a2 a1 a0];

x0_PD1 = [0;0;0]; %Pade 1st Order Initial Conditions
[A1,B1,C1,D1] = tf2ss(b,a); %Transform to State Space
X_PD1( :,1) = x0_PD1; % Initialise with Pade I.Cs
T_PD1(1) = 0; %Time

%2nd Order Coeffs

b4 = 0; b3 = 6;  b2 = 60+18*p;  b1 = 18*p^2+180*p; b0 = 180*p^2;

a4 = 1; a3 = 3*p+8; a2 = 3*p^2+24*p+10;  a1 = 24*p^2-150*p; a0 = 30*p^2;

b = [b4 b3 b2 b1 b0];
a = [a4 a3 a2 a1 a0];
x0_PD2 = [0;0;0;0]; %Pade 2nd Order Initial Conditions
[A2,B2,C2,D2] = tf2ss(b,a); %Transform to State Space
X_PD2( :,1) = x0_PD2; % Initialise with PD I.Cs
T_PD2(1) = 0; % Time

%Backward Euler for Pade Approximation
for n=1:Variable.T/dt
    T_PD1(n+1) = n*dt;
    COEF1 = inv(eye(3)-dt*A1); COEF2 = inv(eye(4)-dt*A2);

    %1st Order Pade
    X_PD1( :,n+1) = COEF1*X_PD1( :,n) + COEF1*dt*B1*Variable.u;

    %2nd Order Pade
    X_PD2( :,n+1) = COEF2*X_PD2( :,n) + COEF2*dt*B2*Variable.u;

end

% Padé Outputs ✅️
Y_PD1 = C1*X_PD1  + D1*Variable.u; Y_PD2 = C2*X_PD2 + D2*Variable.u; Y_OUT_PD = [Y_PD1; Y_PD2];

figure;
plot(T_PD1,[Y; Y_OUT_PD; TD.Y; PD1.Y; PD2.Y] );
grid on;
title('1st and 2nd Order Padé Approximation Delays Vs True Delay');
xlabel('Time (s)');
ylabel('Temperature [°C]');
legend('True Delay','1st Order Pade', '2nd Order Pade', 'True Delay Sim', '1st Order Pade Sim', '2nd Order Pade Sim');
axis ([0 10 0 250])
