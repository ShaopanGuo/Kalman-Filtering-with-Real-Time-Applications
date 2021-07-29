function Xp1 = parrallelAlgEx2(theta0,varTheta0,ite,Q,R,v)

Xp1 = zeros(3,ite+1);
% Xp1(:,1) = [100;100;theta0];
Xp1(:,1) = [0;0;theta0];

P10 = zeros(3,3*(ite+1));
P10(:,1:3) = diag([1,1,varTheta0]);
Gamma = eye(2);

Xp2 = zeros(2,ite+1);
Xp2(:,1) = Xp1(1:2,1);

P20 = zeros(2,2*(ite+1));
P20(:,1:2) = eye(2);

% Q2 = w*0.1*eye(2);
H2 = eye(2);
C2 = [1 0];

for k = 2:ite+1
    % Algorithm 1
    Xkk_1 = [Xp2(1,k-1)+Xp2(2,k-1); Xp1(3,k-1)*Xp2(2,k-1);Xp1(3,k-1)];    % Note this equation
    
    A11 = [1 1; 0 Xp1(3,k-1)];
    A12 = [0; Xp2(2,k-1)];
    A21 = [0 0];
    A22 = 1;
    Ak_1 = [A11 A12; A21 A22];
    GQG = [Gamma*Q*Gamma', zeros(2,1); zeros(1,2), 0];
    Pkk_1 = Ak_1*P10(:,(3*(k-1)-2):3*(k-1))*Ak_1' + GQG;
    
    Gk = Pkk_1*[C2 0]'*inv([C2 0]*Pkk_1*[C2 0]' + R);
    P10(:,(3*k-2):3*k) = (eye(3)-Gk*[C2 0])*Pkk_1;
    Xp1(:,k) = Xkk_1 + Gk*(v(k)-C2*Xkk_1(1:2));
    
    % Algorithm 2
    Fk_1 = [1 1; -0 Xp1(3,k-1)];
    Pkk_1 = Fk_1*P20(:,(2*(k-1)-1):2*(k-1))*Fk_1' + H2*Q*H2';
    Xkk_1 = [Xp2(1,k-1)+Xp2(2,k-1); Xp1(3,k-1)*Xp2(2,k-1)]; 
    Gk = Pkk_1*C2'*inv(C2*Pkk_1*C2' + R);
    P20(:,(2*k-1):2*k) = (eye(2)-Gk*C2)*Pkk_1;
    Xp2(:,k) = Xkk_1 + Gk*(v(k)-C2*Xkk_1);
end
end

