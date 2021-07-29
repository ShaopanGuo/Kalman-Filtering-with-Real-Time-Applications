%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 13.1 The Kalman Smoother 
% Author: Shaopan Guo
% Date: 7/28/2021
% Update: 7/28/2021
% Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Design the filter
% This part is from the official tutorial provided by MATHWORKS
% For this example, set G = B, set H = 0

A = [1.1269   -0.4940    0.1129 
     1.0000         0         0 
          0    1.0000         0];

B = [-0.3832
      0.5919
      0.5191];
  
Gamma = B;

C = [1 0 0];

D = 0;



Ts = -1; % Set the sample time to -1 to mark the plant as discrete 
         % (without a specific sample time).

% Plant dynamics and additive input noise w
sys = ss(A,[B B],C,D,Ts,'InputName',{'u' 'w'},'OutputName','y');  

% The process noise covariance Q and the sensor noise covariance R
Q = 2.3; 
R = 1; 

[kalmf,L,~,Mx,Z] = kalman(sys,Q,R);

sys.InputName = {'u','w'};
sys.OutputName = {'yt'};
vIn = sumblk('y=yt+v');

kalmf.InputName = {'u','y'};
kalmf.OutputName = 'ye';

SimModel = connect(sys,vIn,kalmf,{'u','w','v'},{'yt','ye'});

% generate a known sinusoidal input vector
t = (0:100)';
u = sin(t/5);

% Generate process noise and sensor noise vectors
rng(10,'twister');
w = sqrt(Q)*randn(length(t),1);
v = sqrt(R)*randn(length(t),1);

% simulate the response
out = lsim(SimModel,[u,w,v]);

% generates the response at the outputs yt and ye to the inputs applied 
% at w, v, and u.
yt = out(:,1);   % true response
ye = out(:,2);  % filtered response
y = yt + v;     % measured response


%% Design the Kalman Smoother
%% Based on Ch13.1
N = size(t,1)-1;
% Pf = zeros(size(t,1)*size(A,1),size(A,1));

xs = zeros(size(t,1),size(A,1));
x = zeros(size(t,1),size(A,1));
% xf = x;

for K=1:(N)
    Pf = zeros(size(A));
    Pb = zeros(size(A));
    xf = zeros(size(A,1),1);
    xb = zeros(size(A,1),1);
    for k = 1 : K
%         Pfk_1 = A*Pf(k*size(A,1),size(A,1))*A'+B*Q*B';
%         Gf = Pfk_1*C'/(C*Pfk_1*C'+R);
%         Pf((k*size(A,1)+1):(k+1)*size(A,1),1:size(A,1)) = ...
%             (eye(size(A,1))-Gf*C)*Pfk_1;
        Pfminus1 = A*Pf*A'+B*Q*B';
        Gf = Pfminus1*C'/(C*Pfminus1*C'+R);
        Pf =(eye(size(A,1))-Gf*C)*Pfminus1;
%         xf(k+1,:) = (A*xf(k,:)' + B*u(k) ...
%             + Gf*(y(k+1)-C*A*xf(k,:)' - D*u(k)))';
        xf = A*xf + B*u(k) ...
            + Gf*(y(k+1)-C*A*xf - D*u(k));
    end
    
    for j=1:(N-K)    % j = N-k
        Pbplus1 = Pb + C'/R*C;
        Gb = Pbplus1*Gamma/(Gamma'*Pbplus1*Gamma+inv(Q));
        Pb = A'*(eye(size(A,1))-Gb*Gamma')*Pbplus1*A;
        xb = A'*(eye(size(A,1))-Gb*Gamma')*(xb + C'/R*y(N-j+1)...
            -(C'/R*D + Pbplus1*B)*u(N-j+1));
    end
    
    GK = Pf*Pb/(eye(size(A,1))+Pf*Pb);
    PK = (eye(size(A,1))-GK)*Pf;
    xs(K+1,:) = (eye(size(A,1))-GK)*xf + PK*xb;
end
% Mdl = ssm(A, [B B], C, [D 1]);
% % Mdl = ssm(A, [B B], C, D);
% ys = smooth(Mdl,y);

% Compare the true response with the filtered response.
clf
subplot(211), plot(t,yt,'b',t,ye,'r--',t,xs(:,1),'k--'), 
xlabel('Number of Samples'), ylabel('Output')
title('Kalman Filter Response')
legend('True','Filtered','Smoothed')
subplot(212), plot(t,yt-y,'g',t,yt-ye,'r--',t,yt-xs(:,1),'k--'),
xlabel('Number of Samples'), ylabel('Error')
legend('True - measured','True - filtered','True - smoothed')


