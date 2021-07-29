%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.6 Time-Varing Parameter Identification
% Example 2: System Parameter Estimation
% Author: Shaopan Guo
% Date: 5/13/2021
% Update: 5/20/2021
% Version: 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc

ite = 150;    % Number of iterations

%% Actual Values
X = zeros(3,ite+1);

% Initial values
x0 = 0;
y0 = 0;
z0 = 1.0;        % Here z equals to theta in Ex2 
X(:,1) = [x0; y0; z0];

% Measurements
v = zeros(ite+1,1);

w = 1e-5;
varTheta0 = [0.1, 1.0, 50];
Qt0 = w*diag([0.01,0.01,w]);
% Qt1 = w*diag([0.01,0.01,varTheta(1)]);        % variance of the system noise
% Qt2 = w*diag([0.01,0.01,varTheta(1)]);
% Qt3 = w*diag([0.01,0.01,varTheta(1)]);
Q = w*diag([0.1,0.1]);

H = eye(3);

R = w*0.01;              % variance of the measurement noise

C = [1 0 0];

for k=1:ite
    if k>20 && ~(k>40)
        X(3,k) = 1.2;    
    end
    if k>40 && ~(k>60)
        X(3,k) = 1.5;
    end
    if k>60 && ~(k<80)
        X(3,k) = 2.0;
    end
    if k>80
        X(3,k) = 2.3;
    end
    X(:,k+1) = [1 1 0; 0 X(3,k) 0; 0 0 1]*X(:,k) + sqrtm(Qt0)*randn(3,1);
    v(k) = C*X(:,k) + sqrtm(R)*randn;
end

%% Standard EKF
Xerf = zeros(3,ite+1);
Xerf0 = [100;100;1];
Xerf(:,1) = Xerf0; 

Perf = zeros(3,3*(ite+1));
Perf(:,1:3) = diag([1,1,varTheta0(3)]);

for k=2:ite+1
    Ak_1 = [1 1 0; 0 Xerf(3,k-1) Xerf(2,k-1); 0 0 1];
    Pkk_1 = Ak_1*Perf(:,(3*(k-1)-2):3*(k-1))*Ak_1' + H*Qt0*H';
    Xkk_1 = [Xerf(1,k-1)+Xerf(3,k-1)*Xerf(2,k-1); -0.1*Xerf(1,k-1)+Xerf(2,k-1);Xerf(3,k-1)];
    Gk = Pkk_1*C'*inv(C*Pkk_1*C' + R);
    Perf(:,(3*k-2):3*k) = (eye(3)-Gk*C)*Pkk_1;
    Xerf(:,k) = Xkk_1 + Gk*(v(k)-C*Xkk_1);
end

%% Parallel Alogrithm
theta0 = 5.0;

Xp11 = parrallelAlgEx2(theta0,varTheta0(1),ite,Q,R,v);
Xp12 = parrallelAlgEx2(theta0,varTheta0(2),ite,Q,R,v);
Xp13 = parrallelAlgEx2(theta0,varTheta0(3),ite,Q,R,v);

%% Plot

LW = 2;

figure(811)
kk = 60;
plot(X(3,1:kk),'r-.','LineWidth',LW)
hold on
plot(Xerf(3,1:kk),'b--','LineWidth',LW)
grid on
xlim([0,ite])
legend('actural $$\theta_k$$', 'EKF $$\hat{theta}_k$$', 'Interpreter', 'latex', 'Location','NorthEast')

figure(810)
plot(X(3,:),'r-.','LineWidth',LW)
hold on
plot(Xp11(3,:),'-','LineWidth',LW)
plot(Xp12(3,:),'.','LineWidth',LW)
plot(Xp13(3,:),'--','LineWidth',LW)
grid on
xlim([0,ite])
legend('actual', '$$Var(\theta_0)=0.1$$', '$$Var(\theta_0)=1.0$$', '$$Var(\theta_0)=50$$', 'Interpreter', 'latex', 'Location','NorthEast')