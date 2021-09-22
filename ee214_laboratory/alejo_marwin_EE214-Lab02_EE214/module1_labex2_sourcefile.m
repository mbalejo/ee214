%% Marwin B. Alejo   2020-20221   EE214_Module1-LabEx2
% * Date Performed (d/m/y): 22/09/2021

%% I. Binary Symmetric Channel
% 1. (Manual Computation) Consider a random transmission of 100 bits on a BSC with a 1-p=0.05
% and p=0.95, determine the probabiities of the ff. Let T = Transmitted and
% R = Received.

N_100 = 100; % number of bits to be generated
tx_100 = rand(1,N_100) > 0.5; % generate 100 random bits
tx0_100 = (sum(tx_100(:)==0))/100; % number of transmitted 0's
tx1_100 = (sum(tx_100(:)==1))/100; % number of transmitted 1's
txe_100 = 0.05; % transmission error probability rate
txs_100 = 0.95; % transmission success probability rate 

%% 1. Manual Computation of P()
% To determine error and success P(Rx) and P(Tx) through MANUAL 
% computations, the ff. equations must be considered:

%% a.1 P(R'1') receive P() success rate
% 
% $$P(R'1')=P(T'1')\times P(R'1'\mid T'1')+P(T'0')\times P(R'1'\mid T'0')$$
% 
rx1_100 = tx1_100*txs_100+tx0_100*txe_100; % computes the P(R'1')
fprintf('P(R1) is equal to %.4f',rx1_100); % display P(R'1')

% a.1. To determine P(R'1') and P(R'0') through MANUAL computations, the ff. 
% equations must be considered:

%% a.2 P(R'0') receive P() error rate
% 
% $$P(R'0')=P(T'0')\times P(R'0'\mid T'0')+P(T'1')\times P(R'0'\mid T'1')$$
% 
rx0_100 = tx0_100*txs_100+tx1_100*txe_100; % computes the P(R'0')
fprintf('P(R0) is equal to %.4f',rx0_100); % display P(R'0')

%% b.1 P(R'1'|T'1') - conditional P() aka priti
% 
% $$P(R'1'|T'1')=\frac{P(T'1'|R'1')\times P(R'1')}{P(T'1')}$$
% 
PR1T1_100 = (txs_100*rx1_100)/tx1_100; % computes P(R'1'|T'1')
fprintf('P(R1|T1) is equal to %.4f', PR1T1_100); %display P(R'1'|T'1')

%% b.2 P(R'0'|T'1') - conditional P() aka proti
% 
% $$P(R'0'|T'1')=\frac{P(T'1'|R'0')\times P(R'0')}{P(T'1')}$$
% 
PR0T1_100 = (txe_100*rx0_100)/tx1_100; % computes P(R'0'|T'1')
fprintf('P(R0|T1) is equal to %.4f', PR0T1_100); %display P(R'0'|T'1')

%% b.3 P(R'1'|T'0') - conditional P() aka prito
% 
% $$P(R'1'|T'0')=\frac{P(T'0'|R'1')\times P(R'1')}{P(T'0')}$$
% 
PR1T0_100 = (txe_100*rx1_100)/tx0_100; % computes P(R'1'|T'0')
fprintf('P(R1|T0) is equal to %.4f', PR1T0_100); %display P(R'1'|T'0')

%% b.4 P(R'0'|T'0') - conditional P() aka proto
% 
% $$P(R'0'|T'0')=\frac{P(T'0'|R'0')\times P(R'0')}{P(T'0')}$$
% 
PR0T0_100 = (txs_100*rx0_100)/tx0_100; % computes P(R'0'|T'0')
fprintf('P(R0|T0) is equal to %.4f', PR0T0_100); %display P(R'0'|T'0')

%% c.1 P(T'1'|R'1') - bayesian P() aka P'tiri
% 
% $$P(T'1'|R'1')=\frac{P(R'1'|T'1')\times P(T'1')}{P(R'1'|T'1')\times P(T'1') + P(R'1'|T'0')\times P(T'0')}$$
% 
PT1R1_100 = (txs_100*tx1_100)/((txs_100*tx1_100)+(txe_100*tx0_100)); % computes P(T'1'|R'1')
fprintf('P(T1|R1) is equal to %.4f', PT1R1_100); %display P(T'1'|R'1')

%% c.2 P(T'0'|R'1') - bayesian P() aka P'tori
% 
% $$P(T'0'|R'1')=\frac{P(R'1'|T'0')\times P(T'0')}{P(R'1'|T'0')\times P(T'0') + P(R'1'|T'1')\times P(T'1')}$$
% 
PT0R1_100 = (txe_100*tx0_100)/((txe_100*tx0_100)+(txs_100*tx1_100)); % computes P(T'0'|R'1')
fprintf('P(T0|R1) is equal to %.4f', PT0R1_100); %display P(T'0'|R'1')

%% c.3 P(T'1'|R'0') - bayesian P() aka P'tiro
% 
% $$P(T'1'|R'0')=\frac{P(R'0'|T'1')\times P(T'1')}{P(R'0'|T'1')\times P(T'1') + P(R'0'|T'0')\times P(T'0')}$$
% 
PT1R0_100 = (txe_100*tx1_100)/((txe_100*tx1_100)+(txs_100*tx0_100)); % computes P(T'1'|R'0')
fprintf('P(T1|R0) is equal to %.4f', PT1R0_100); %display P(T'1'|R'0')

%% c.4 P(T'0'|R'0') - bayesian P() aka P'toro
% 
% $$P(T'1'|R'0')=\frac{P(R'0'|T'1')\times P(T'1')}{P(R'0'|T'1')\times P(T'1') + P(R'0'|T'0')\times P(T'0')}$$
% 
PT1R0_100 = (txe_100*tx1_100)/((txe_100*tx1_100)+(txs_100*tx0_100)); % computes P(T'1'|R'0')
fprintf('P(T1|R0) is equal to %.4f', PT1R0_100); %display P(T'1'|R'0')
