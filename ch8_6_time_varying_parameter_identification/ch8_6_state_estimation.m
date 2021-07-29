%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.6 Time-Varing Parameter Identification
% Example 1: State Estimation
% Author: Shaopan Guo
% Date: 5/13/2021
% Update: 5/20/2021
% Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc

ite = 500;    % Number of iterations

%% Uncorrelated zero-mean Gaussian white noise
% zeta and eta are uncorrelated zero-mean Gaussian white noise sequences
% with Var(zeta(k))=0.1*I3 and Var(eta(k))=0.01   

%% Actual Values
X = zeros(3,ite+1);

% Initial values
x0 = 1.0;
y0 = 1.0;
z0 = 0.1;
X(:,1) = [x0; y0; z0];

% Measurements
v = zeros(ite+1,1);

w = 1e-5;
Q = w*0.1*eye(3);        % variance of the system noise

H = eye(3);

R = w*0.01;              % variance of the measurement noise


C = [1 0 0];

for k=1:ite
    X(:,k+1) = [1 X(3,k) 0; -0.1 1 0; 0 0 1]*X(:,k) + sqrtm(Q)*randn(3,1);
    v(k) = C*X(:,k) + sqrtm(R)*randn;
end

%% Standard EKF
Xerf = zeros(3,ite+1);
Xerf0 = [100;100;1];
Xerf(:,1) = Xerf0; 

Perf = zeros(3,3*(ite+1));
Perf(:,1:3) = eye(3);

for k=2:ite+1
    Ak_1 = [1 Xerf(3,k-1) Xerf(2,k-1); -0.1 1 0; 0 0 1];
    Pkk_1 = Ak_1*Perf(:,(3*(k-1)-2):3*(k-1))*Ak_1' + H*Q*H';
    Xkk_1 = [Xerf(1,k-1)+Xerf(3,k-1)*Xerf(2,k-1); -0.1*Xerf(1,k-1)+Xerf(2,k-1);Xerf(3,k-1)];
    Gk = Pkk_1*C'*inv(C*Pkk_1*C' + R);
    Perf(:,(3*k-2):3*k) = (eye(3)-Gk*C)*Pkk_1;
    Xerf(:,k) = Xkk_1 + Gk*(v(k)-C*Xkk_1);
end

%% Parallel Alogrithm

Xp1 = zeros(3,ite+1);
Xp1(:,1) = Xerf0; 

P10 = zeros(3,3*(ite+1));
P10(:,1:3) = eye(3);

Xp2 = zeros(2,ite+1);
Xp2(:,1) = Xerf0(1:2); 

P20 = zeros(2,2*(ite+1));
P20(:,1:2) = eye(2);

Q2 = w*0.1*eye(2);
H2 = eye(2);
C2 = [1 0];

for k = 2:ite+1
    % Algorithm 1
    Ak_1 = [1 Xp1(3,k-1) Xp2(2,k-1); -0.1 1 0; 0 0 1];
    Pkk_1 = Ak_1*P10(:,(3*(k-1)-2):3*(k-1))*Ak_1' + H*Q*H';
    Xkk_1 = [Xp2(1,k-1)+Xp1(3,k-1)*Xp2(2,k-1); -0.1*Xp2(1,k-1)+Xp2(2,k-1);Xp1(3,k-1)];    % Note this equation
    Gk = Pkk_1*C'*inv(C*Pkk_1*C' + R);
    P10(:,(3*k-2):3*k) = (eye(3)-Gk*C)*Pkk_1;
    Xp1(:,k) = Xkk_1 + Gk*(v(k)-C*Xkk_1);
    
    % Algorithm 2
    Fk_1 = [1 Xp1(3,k-1); -0.1 1];
    Pkk_1 = Fk_1*P20(:,(2*(k-1)-1):2*(k-1))*Fk_1' + H2*Q2*H2';
    Xkk_1 = [Xp2(1,k-1)+Xp1(3,k-1)*Xp2(2,k-1); -0.1*Xp2(1,k-1)+Xp2(2,k-1)]; 
    Gk = Pkk_1*C2'*inv(C2*Pkk_1*C2' + R);
    P20(:,(2*k-1):2*k) = (eye(2)-Gk*C2)*Pkk_1;
    Xp2(:,k) = Xkk_1 + Gk*(v(k)-C2*Xkk_1);
end


%% Plot

LW = 2;

figure(84)
plot(X(1,:),'r-.','LineWidth',LW)
hold on
plot(Xerf(1,:),'b--','LineWidth',LW)
grid on
xlim([0,ite])
legend('actural $$x_k$$', 'EKF $$\hat{x}_k$$', 'Interpreter', 'latex', 'Location','SouthWest')

figure(85)
plot(X(1,:),'r-.','LineWidth',LW)
hold on
plot(Xp1(1,:),'k--','LineWidth',LW)
grid on
xlim([0,ite])
legend('actural $$x_k$$', 'parallel algorithm $$\tilde{x}_k$$', 'Interpreter', 'latex', 'Location','SouthWest')

figure(86)
plot(X(2,:),'r-.','LineWidth',LW)
hold on
plot(Xerf(2,:),'b--','LineWidth',LW)
grid on
xlim([0,ite])
legend('actural $$y_k$$', 'EKF $$\hat{y}_k$$', 'Interpreter', 'latex', 'Location','SouthWest')

figure(87)
plot(X(2,:),'r-.','LineWidth',LW)
hold on
plot(Xp1(2,:),'k--','LineWidth',LW)
grid on
xlim([0,ite])
legend('actural $$y_k$$', 'parallel algorithm $$\tilde{y}_k$$', 'Interpreter', 'latex', 'Location','SouthWest')

figure(88)
plot(X(3,:),'r-.','LineWidth',LW)
hold on
plot(Xerf(3,:),'b--','LineWidth',LW)
grid on
xlim([0,ite])
legend('actural $$z_k$$', 'EKF $$\hat{z}_k$$', 'Interpreter', 'latex', 'Location','SouthWest')

figure(89)
plot(X(3,:),'r-.','LineWidth',LW)
hold on
plot(Xp1(3,:),'k--','LineWidth',LW)
grid on
xlim([0,ite])
legend('actural $$z_k$$', 'parallel algorithm $$\tilde{z}_k$$', 'Interpreter', 'latex', 'Location','SouthWest')

figure(90)
plot(1:size(P10(1,1:2:3*ite),2),P10(1,1:2:3*ite),'r-.','LineWidth',LW)
hold on
plot(1:size(Perf(1,1:2:3*ite),2),Perf(1,1:2:3*ite),'b-.','LineWidth',LW)
grid on
legend('parallel algorithm', 'EKF')

figure(91)
plot(1:size(P10(2,2:2:3*ite),2),P10(2,2:2:3*ite),'r-.','LineWidth',LW)
hold on
plot(1:size(Perf(2,2:2:3*ite),2),Perf(2,2:2:3*ite),'b-.','LineWidth',LW)
grid on
legend('parallel algorithm', 'EKF')

figure(92)
plot(1:size(P10(3,3:2:3*ite),2),P10(3,3:2:3*ite),'r-.','LineWidth',LW)
hold on
plot(1:size(Perf(3,3:2:3*ite),2),Perf(3,3:2:3*ite),'b-.','LineWidth',LW)
grid on
legend('parallel algorithm', 'EKF')

figure(93)
plot(1:size(P10(1,1:2:3*ite),2),P10(1,1:2:3*ite)-Perf(1,1:2:3*ite),'r-.','LineWidth',LW)
hold on
grid on

figure(94)
plot(1:size(P10(2,2:2:3*ite),2),P10(2,2:2:3*ite)-Perf(2,2:2:3*ite),'r-.','LineWidth',LW)
hold on
grid on

figure(95)
plot(1:size(P10(3,3:2:3*ite),2),P10(3,3:2:3*ite)-Perf(3,3:2:3*ite),'r-.','LineWidth',LW)
hold on
grid on

P30 = P10-Perf;
for k=1:ite
    eig(P30(:,3*(k-1)+1:3*k))
end