%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9.2 Real-Time Tracking
% Equation 9.12.3-9.12.5: Verify the Cramer's Rule used in Hi 
% Author: Shaopan Guo
% Date: 6/9/2021
% Update: 6/10/2021
% Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms z h v g1 g2 g3

Delta = [z-1+g1, -h+h*g1, -h^2/2+h^2*g1/2;
    g2, z-1+h*g2, -h+h^2*g2/2;
    g3, h*g3, z-1+h^2*g3/2];
det(Delta);

B = z * [g1; g2; g3] * v;

Delta1 = [B, Delta(:,2:3)];
Delta2 = [Delta(:,1), B, Delta(:,3)];
Delta3 = [Delta(:,1:2), B];

H1 = det(Delta1)/(v*det(Delta));
H2 = det(Delta2)/(v*det(Delta));
H3 = det(Delta3)/(v*det(Delta));

% pretty(H1)
% pretty(H2)
% pretty(H3)

% Final results of the book are right
pretty(det(Delta))
pretty(det(Delta1))
pretty(det(Delta2))
pretty(det(Delta3))

