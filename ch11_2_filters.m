%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11.1.2 Discrete Wavelet and Filter Banks
% Calculate filter operators and verify their properties
% Author: Shaopan Guo
% Date: 7/19/2021
% Update: 7/21/2021
% Version: 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variable List
% H1: H^{i-1}
% H2: H^{i-2}
% G1: G^{i-1}
% G2: G^{i-2}

%% Calculate H^{i-1} and G^{i-1} for M=4

lowPassFilter4 = sqrt(2)/2 * [1 1 0 0;
                              0 1 1 0;
                              0 0 1 1;
                              0 0 0 1];
highPassFilter4 = sqrt(2)/2 * [ 1 -1  0  0;
                                0  1 -1  0;
                                0  0  1 -1;
                                0  0  0  1];
downSampling4 = [1 0 0 0;
                 0 0 1 0];
          
H1 = downSampling4*lowPassFilter4;
G1 = downSampling4*highPassFilter4;

%% Verify the properties of the filters

H1' * H1 + G1' * G1

[H1*H1' H1*G1' ; G1*H1' G1*G1']


%% Calculate H^{i-2} and G^{i-2} for M=4

lowPassFilter2 = sqrt(2)/2 * [1 1; 
                              0 1];
highPassFilter2 = sqrt(2)/2 * [ 1 -1; 
                                0  1];
downSampling2 = [1 0];
          
H2 = downSampling2*lowPassFilter2;
G2 = downSampling2*highPassFilter2;

%% Verify the properties of the filters

H2' * H2 + G2' * G2

[H2*H2' H2*G2' ; G2*H2' G2*G2']
            
