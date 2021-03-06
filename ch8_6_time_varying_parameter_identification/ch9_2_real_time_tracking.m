%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9.2 Real-Time Tracking 
% Author: Xinxin Wang
% Reviewer: Shaopan Guo
% Date: 6/10/2021
% Update: 6/10/2021
% Version: 1.0
% Report: No bugs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter Table

% h: sampling time 0.2
% x: real state values x_k \dot{x}_k \ddot{x}_k
% x_e: estimated state values
% x1; disturbed state values
% z0: measurement

%%
clc; clear; close all;

%% init values
h = 0.2;               
x(:, 3)=[4; 5; 6];      % x_{0} = x(:,3)
% x_e(:, 3)=[4; 5; 6];    %                        No initial estimate error?
x_e(:, 3)=[0; 0; 0];    % It does not matter.
                        %                    What does following part mean?
                        % The following part gives x_k for k<-1
                        %                    x_k = 0 according to the book?
for k=1:2
%     x_e(:, k)=[1;1;1];    % It does not matter.
    x_e(:, k)=[0;0;0];
%     z0(k)=1;            %                                What does z0 mean?
                        % z0 is the measurement. 
                        % It should be 0 according to the book.
                        % It does not matter.
    z0(k)=0;    
end

%% System
% See the help document for kalman to get the definitions of A,B,C,D,G,H
A=[1 h h^2/2;0 1 h;0 0 1];      
B=zeros(3, 3);
G=eye(3, 3);                    %                         What does G mean?
                                % G is the input matrix for the sys. noise
C=[1, 0, 0];
D=[0, 0, 0];
H=zeros(1, 3);                  %                         What does H mean?
                                % H is the input matrix for the obs. noise
                                
for k = 4:50
    % model noise ~ N(0, 1)
    w=randn(3, 1);
    % obs noise ~ N(0, 1)
    v=randn(1);
    x(:, k)=A*x(:, k-1);
    % add model noise
    x1(:, k)=x(:, k)+w;          
    QN=eye(3, 3);
    RN=1;
    NN=0;
    % add obs noise
    z0(k)=C*x1(:, k)+v;
   
    % SYS=ss(A, [B G], C, [D H], 1);                    % Here 1 should be h?
    SYS=ss(A, [B G], C, [D H], h);                   
    [KEST, L, P]=kalman(SYS, QN, RN, NN);
    % get filteration
    %x_e(:, k)=A*x_e(:, k-1)+L.*(z0(:, k)-C*A*x_e(:, k-1));
    x_e(1, k)=-((L(1,1)-3)+L(2,1)*h+L(3,1)*h^2/2)*x_e(1, k-1)-((3-2*L(1,1))-L(2,1)*h+L(3,1)*h^2/2)*x_e(1, k-2)-(L(1,1)-1)*x_e(1, k-3)+L(1,1)*z0(k)+(L(3,1)*h^2/2+L(2,1)*h-2*L(1,1))*z0(k-1)+(L(3,1)*h^2/2-L(2,1)*h+L(1,1))*z0(k-2);
    x_e(2, k)=-((L(1,1)-3)+L(2,1)*h+L(3,1)*h^2/2)*x_e(2, k-1)-((3-2*L(1,1))-L(2,1)*h+L(3,1)*h^2/2)*x_e(2, k-2)-(L(1,1)-1)*x_e(2, k-3)+L(2,1)*z0(k)+(h*L(3,1)-2*L(2,1))*z0(k-1)+(L(2,1)-h*L(3,1))*z0(k-2);
    x_e(3, k)=-((L(1,1)-3)+L(2,1)*h+L(3,1)*h^2/2)*x_e(3, k-1)-((3-2*L(1,1))-L(2,1)*h+L(3,1)*h^2/2)*x_e(3, k-2)-(L(1,1)-1)*x_e(3, k-3)+L(3,1)*z0(k)-2*L(3,1)*z0(k-1)+L(3,1)*z0(k-2); 

end

%% plot and compare
k=1:50;
plot(k,x(1,k),'r');
hold on
plot(k,x_e(1,k),'r-.');
hold on
plot(k,x(2,k),'b');
hold on
plot(k,x_e(2,k),'b-.');
hold on
plot(k,x(3,k),'g');
hold on
plot(k,x_e(3,k),'g-.');
legend('real value of x[1]','estimate of x[1]',  'real value of x[2]','estimate of x[2]', 'real value of x[3]','estimate of x[3]');
xlabel('k');
ylabel('x');
