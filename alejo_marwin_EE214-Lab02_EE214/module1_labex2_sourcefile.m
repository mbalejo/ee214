%% Marwin B. Alejo   2020-20221   EE214_Module1-LabEx2
% * Date Performed (d/m/y): 22/09/2021
% * Disclaimer: values of the used variables varies for every execurtion as
% these are dependent on the randomly defined values of 1s and 0s in N, 
% answers that require particular values are directed to the
% variable where it is stored.

%% I. Binary Symmetric Channel
% Consider a random transmission of 100 bits on a BSC with a 1-p=0.05
% and p=0.95, determine the probabiities of the ff. Let T = Transmitted and
% R = Received.

N_100 = 100; % number of bits to be generated
tx_100 = rand(1,N_100) > 0.5; % generate 100 random bits
tx0_100 = (sum(tx_100(:)==0))/N_100 % rate of transmitted 0's from 100 sample bits
tx1_100 = (sum(tx_100(:)==1))/N_100 % rate of transmitted 1's from 100 sample bits
txe_100 = 0.05; % transmission error probability rate
txs_100 = 0.95; % transmission success probability rate 

%% 1. Manual Computation of P()
% To determine error and success P(Rx) and P(Tx) through MANUAL 
% computations, the ff. equations must be considered:

% Manually computed P(R) - vector | manual XOR implementation
manual_rx_100 = tx_100.*txs_100+tx_100.*txe_100;

 %% a.1 P(R'1'): P() rxing 1s
% 
% $$P(R'1')=P(T'1')\times P(R'1'\mid T'1')+P(T'0')\times P(R'1'\mid T'0')$$
% 
rx1_100 = tx1_100*txs_100+tx0_100*txe_100; % computes the P(R'1')
fprintf('Given the conditions above, the probability that the receiver receive 1s from 100 sample bits is %.4f',rx1_100); % display P(R'1')

%% a.2 P(R'0'): P() rxing 0s
% 
% $$P(R'0')=P(T'0')\times P(R'0'\mid T'0')+P(T'1')\times P(R'0'\mid T'1')$$
% 
rx0_100 = tx0_100*txs_100+tx1_100*txe_100; % computes the P(R'0')
fprintf('Given the conditions above, the probability that the receiver receive 0s from 100 sample bits is %.4f',rx0_100); % display P(R'0')

%% b.1 P(R'1'|T'1') - conditional P() aka priti
% 
% $$P(R'1'|T'1')=\frac{P(T'1'|R'1')\times P(R'1')}{P(T'1')}$$
% 
PR1T1_100 = (txs_100*rx1_100)/tx1_100; % computes P(R'1'|T'1')
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 1 is %.4f', PR1T1_100); %display P(R'1'|T'1')

%% b.2 P(R'0'|T'1') - conditional P() aka proti
% 
% $$P(R'0'|T'1')=\frac{P(T'1'|R'0')\times P(R'0')}{P(T'1')}$$
% 
PR0T1_100 = (txe_100*rx0_100)/tx1_100; % computes P(R'0'|T'1')
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitter transmitted bit 1 is %.4f', PR0T1_100); %display P(R'0'|T'1')

%% b.3 P(R'1'|T'0') - conditional P() aka prito
% 
% $$P(R'1'|T'0')=\frac{P(T'0'|R'1')\times P(R'1')}{P(T'0')}$$
% 
PR1T0_100 = (txe_100*rx1_100)/tx0_100; % computes P(R'1'|T'0')
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 0 is %.4f', PR1T0_100); %display P(R'1'|T'0')

%% b.4 P(R'0'|T'0') - conditional P() aka proto
% 
% $$P(R'0'|T'0')=\frac{P(T'0'|R'0')\times P(R'0')}{P(T'0')}$$
% 
PR0T0_100 = (txs_100*rx0_100)/tx0_100; % computes P(R'0'|T'0')
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitter transmitted bit 0 is %.4f', PR0T0_100); %display P(R'0'|T'0')

%% c.1 P(T'1'|R'1') - bayesian P() aka P'tiri
% 
% $$P(T'1'|R'1')=\frac{P(R'1'|T'1')\times P(T'1')}{P(R'1'|T'1')\times P(T'1') + P(R'1'|T'0')\times P(T'0')}$$
% 
PT1R1_100 = (txs_100*tx1_100)/((txs_100*tx1_100)+(txe_100*tx0_100)); % computes P(T'1'|R'1')
fprintf('Given the conditions above, the probability that the transmitter transmitted bit 1 when the received data by the receiver is bit 1 is %.4f', PT1R1_100); %display P(T'1'|R'1')

%% c.2 P(T'0'|R'1') - bayesian P() aka P'tori
% 
% $$P(T'0'|R'1')=\frac{P(R'1'|T'0')\times P(T'0')}{P(R'1'|T'0')\times P(T'0') + P(R'1'|T'1')\times P(T'1')}$$
% 
PT0R1_100 = (txe_100*tx0_100)/((txe_100*tx0_100)+(txs_100*tx1_100)); % computes P(T'0'|R'1')
fprintf('Given the conditions above, the probability that the transmitter transmitted bit 0 when then receiver received bit 1 is %.4f', PT0R1_100); %display P(T'0'|R'1')

%% c.3 P(T'1'|R'0') - bayesian P() aka P'tiro
% 
% $$P(T'1'|R'0')=\frac{P(R'0'|T'1')\times P(T'1')}{P(R'0'|T'1')\times P(T'1') + P(R'0'|T'0')\times P(T'0')}$$
% 
PT1R0_100 = (txe_100*tx1_100)/((txe_100*tx1_100)+(txs_100*tx0_100)); % computes P(T'1'|R'0')
fprintf('Given the conditions above, the probability that the transmitter transmitted bit 1 when the received data by the receiver is bit 0 is %.4f', PT1R0_100); %display P(T'1'|R'0')

%% c.4 P(T'0'|R'0') - bayesian P() aka P'toro
% 
% $$P(T'0'|R'0')=\frac{P(R'0'|T'0')\times P(T'0')}{P(R'0'|T'0')\times P(T'0') + P(R'0'|T'1')\times P(T'1')}$$
% 
PT0R0_100 = (txs_100*tx0_100)/((txs_100*tx0_100)+(txe_100*tx1_100)); % computes P(T'0'|R'0')
fprintf('Given the conditions above, the probability that the transmitter send bit 0 when the receiver received bit 0 is %.4f', PT0R0_100); %display P(T'0'|R'0')

%% 2. MATLAB Computation of P()
% To determine error and success P(Rx) and P(Tx) through MATLAB, the ff. 
% codes must be considered:

% For comparison, this figure show the Tx and Rx of both the manually 
% computed (see line 25) and MATLAB-generated Rx vectors.

% lines 13-14 must be inserted here!!!
figure(); subplot(2,2,2)
stem(tx_100);
xlabel('Transmitted bits'); title('MATLAB Computation'); 
% determine is successfully transmitted or not
ch_100 = rand(1,N_100) > txs_100;
rx_100 = xor(tx_100, ch_100);
%plot the results'
subplot(2,2,4), stem(rx_100);
xlabel('Received bits');
subplot(2,2,1), stem(tx_100); xlabel('Transmitted bits'); title('Manual Computation');
subplot(2,2,3), stem(manual_rx_100); xlabel('Received bits'); 
%% 

% P(R'1') and P(R'0') - MATLAB values

% rate of received 0's from 100 sample bits
matlab_rx0_100 = (sum(rx_100(:)==0))/N_100

% rate of received 1's from 100 sample bits
matlab_rx1_100 = (sum(rx_100(:)==1))/N_100
%% 

% PR1T1 from MATLAB gen values
matlab_PR1T1_100 = (txs_100*matlab_rx1_100)/tx1_100 
%% 

% PR0T1 from MATLAB gen values
matlab_PR0T1_100 = (txe_100*matlab_rx0_100)/tx1_100 
%% 

% PR0T1 from MATLAB gen values
matlab_PR1T0_100 = (txe_100*matlab_rx1_100)/tx0_100 
%% 

% PR0T0 from MATLAB gen values
matlab_PR0T0_100 = (txs_100*matlab_rx0_100)/tx0_100 
%% 

% PT1R1 from MATLAB gen values
matlab_PT1R1_100 = (txs_100*tx1_100)/((txs_100*tx1_100)+(txe_100*tx0_100)) 
%% 

% PT0R1 from MATLAB gen values
matlab_PT0R1_100 = (txe_100*tx0_100)/((txe_100*tx0_100)+(txs_100*tx1_100)) 
%% 

% PT1R0 from MATLAB gen values
matlab_PT1R0_100 = (txe_100*tx1_100)/((txe_100*tx1_100)+(txs_100*tx0_100))
%% 

% PT0R0 from MATLAB gen values
matlab_PT0R0_100 = (txs_100*tx0_100)/((txs_100*tx0_100)+(txe_100*tx1_100))

%% Discussion for 1 and 2
% * The calculated P(R'1') and P(R'0') values are defined by the
% variables 'rx1_100' and 'rx0_100' for manual computation (see sec.a1-a2)
% and 'matlab_rx0_100' and 'matlab_rx1_100'. Although the small variance
% between the manually computed and MATLAB-generated received bits are not visible from the figure above, the
% numbers of these variables show that the there exist an error rate of up
% to ~5% (depending on the computed values in line 14) within the 
% transmission channel.
%
% * Similarly, the manually computed (sec.b1-c4) and MATLAB-generated (lines 116-157) 
% PR0T0, PR0T1, PR1T0, PR1T1, T1R1, R1T1, R1T0, and R0T0 show a tiny variance to none from each other assuring
% that there is only a chance of up to ~5% that the received or transmitted bits are
% erroneous. Moreover, the success transmission/receving rate of correct 1
% or 0 bits ccording to these probabilities is up to ~95%, which is aligned
% to the provided specification above.
%
% * By inspection, the values from the manually computed probabilities of the samples by
% conditional differ when the same means applied onto the MATLAB generated
% bits. In contrast to when these probabilities are computed through
% Bayesian theorem, the computed values do not change.

%% 3. Repeat experiment for N=1000 and N=10000 bits.
% To determine error and success P(Rx) and P(Tx) when N=1k and N=10k bits,
% the ff. codes must be realized.

%% 3.a.1 N=1000 bits (Manual Computation)
N_1k = 1000; % number of bits to be generated
tx_1k = rand(1,N_1k) > 0.5; % generate 1000 random bits
tx0_1k = (sum(tx_1k(:)==0))/N_1k % rate of tx 0's from 1k sample bits
tx1_1k = (sum(tx_1k(:)==1))/N_1k % rate of tx 1's from 1k sample bits
txe_1k = 0.05; % transmission error probability rate
txs_1k = 0.95; % transmission success probability rate 

% Manually computed P(R) - vector | manual XOR implementation
manual_rx_1k = tx_1k.*txs_1k+tx_1k.*txe_1k;
%% 

% P(R'1') of 1k bits
rx1_1k = tx1_1k*txs_1k+tx0_1k*txe_1k;
fprintf('Given the conditions above, the probability that the receiver receive 1s from 1000 sample bits is %.4f',rx1_1k);
%% 

% P(R'0') of 1k bits
rx0_1k = tx0_1k*txs_1k+tx1_1k*txe_1k;
fprintf('Given the conditions above, the probability that the receiver receive 0s from 1000 sample bits is %.4f',rx0_1k);
%% 

% P(R'1'|T'1') of 1k bits
PR1T1_1k = (txs_1k*rx1_1k)/tx1_1k;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 1 is %.4f', PR1T1_1k);
%% 

% P(R'0'|T'1') of 1k bits
PR0T1_1k = (txe_1k*rx0_1k)/tx1_1k;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitter transmitted bit 1 is %.4f', PR0T1_1k);
%% 

% P(R'1'|T'0') of 1k bits
PR1T0_1k = (txe_1k*rx1_1k)/tx0_1k;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 0 is %.4f', PR1T0_1k);
%% 

% P(R'0'|T'0') of 1k bits
PR0T0_1k = (txs_1k*rx0_1k)/tx0_1k;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitter transmitted bit 0 is %.4f', PR0T0_1k);
%% 

% P(T'1'|R'1') of 1k bits
PT1R1_1k = (txs_1k*tx1_1k)/((txs_1k*tx1_1k)+(txe_1k*tx0_1k));
fprintf('Given the conditions above, the probability that the transmitter transmitted bit 1 when the received data by the receiver is bit 1 is %.4f', PT1R1_1k);
%% 

% P(T'0'|R'1') of 1k bits
PT0R1_1k = (txe_1k*tx0_1k)/((txe_1k*tx0_1k)+(txs_1k*tx1_1k));
fprintf('Given the conditions above, the probability that the transmitter transmitted bit 0 when then receiver received bit 1 is %.4f', PT0R1_1k);
%% 

% P(T'1'|R'0') of 1k bits
PT1R0_1k = (txe_1k*tx1_1k)/((txe_1k*tx1_1k)+(txs_1k*tx0_1k));
fprintf('Given the conditions above, the probability that the transmitter transmitted bit 1 when the received data by the receiver is bit 0 is %.4f', PT1R0_1k);
%% 

% P(T'0'|R'0') of 1k bits
PT0R0_1k = (txs_1k*tx0_1k)/((txs_1k*tx0_1k)+(txe_1k*tx1_1k));
fprintf('Given the conditions above, the probability that the transmitter send bit 0 when the receiver received bit 0 is %.4f', PT0R0_1k);

%% 3.a.2 N=1000 bits (MATLAB Computation)
%% 

% lines 184-185 must be inserted here!!!
figure(); subplot(4,1,3)
stem(tx_1k);
xlabel('Transmitted 1k bits'); title('MATLAB Computation'); 
% determine is successfully transmitted or not
ch_1k = rand(1,N_1k) > txs_1k;
rx_1k = xor(tx_1k, ch_1k);
%plot the results'
subplot(4,1,4), stem(rx_1k);
xlabel('Received 1k bits');
subplot(4,1,1), stem(tx_1k); xlabel('Transmitted 1k bits'); title('Manual Computation');
subplot(4,1,2), stem(manual_rx_1k); xlabel('Received 1k bits'); 

%% 

% P(R'1') and P(R'0') - MATLAB values

% rate of received 0's from 1k sample bits
matlab_rx0_1k = (sum(rx_1k(:)==0))/N_1k

% rate of received 1's from 1k sample bits
matlab_rx1_1k = (sum(rx_1k(:)==1))/N_1k
%% 

% PR1T1 from MATLAB gen values - 1k bits
matlab_PR1T1_1k = (txs_1k*matlab_rx1_1k)/tx1_1k 
%% 

% PR0T1 from MATLAB gen values - 1k bits
matlab_PR0T1_1k = (txe_1k*matlab_rx0_1k)/tx1_1k 
%% 

% PR0T1 from MATLAB gen values - 1k bits
matlab_PR1T0_1k = (txe_1k*matlab_rx1_1k)/tx0_1k 
%% 

% PR0T0 from MATLAB gen values - 1k bits
matlab_PR0T0_1k = (txs_1k*matlab_rx0_1k)/tx0_1k 
%% 

% PT1R1 from MATLAB gen values - 1k bits
matlab_PT1R1_1k = (txs_1k*tx1_1k)/((txs_1k*tx1_1k)+(txe_1k*tx0_1k)) 
%% 

% PT0R1 from MATLAB gen values - 1k bits
matlab_PT0R1_1k = (txe_1k*tx0_1k)/((txe_1k*tx0_1k)+(txs_1k*tx1_1k)) 
%% 

% PT1R0 from MATLAB gen values - 1k bits
matlab_PT1R0_1k = (txe_1k*tx1_1k)/((txe_1k*tx1_1k)+(txs_1k*tx0_1k))
%% 

% PT0R0 from MATLAB gen values - 1k bits
matlab_PT0R0_1k = (txs_1k*tx0_1k)/((txs_1k*tx0_1k)+(txe_1k*tx1_1k))

%% Discussion for 3a1-2
% * The calculated P(R'1') and P(R'0') values are defined by the
% variables 'rx1_1k' and 'rx0_1k' for manual computation (see sec.3a1)
% and 'matlab_rx0_1k' and 'matlab_rx1_1k'. Although the small variance
% between the manually computed and MATLAB-generated received bits are not visible from the figure above, the
% numbers of these variables show that the there exist an error rate of up
% to ~5% (depending on the computed values of N) within the 
% transmission channel.
%
% * Similarly, the manually computed (sec.b1-c4) and MATLAB-generated (lines 271-300) 
% PR0T0, PR0T1, PR1T0, PR1T1, T1R1, R1T1, R1T0, and R0T0 show a tiny variance to none from each other assuring
% that there is only a chance of up to ~5% that the received or transmitted bits are
% erroneous. Moreover, the success transmission/receving rate of correct 1
% or 0 bits ccording to these probabilities is up to ~95%, which is aligned
% to the provided specification above.
%
% * By inspection, the values from the manually computed probabilities of the samples by
% conditional differ when the same means applied onto the MATLAB generated
% bits. In contrast to when these probabilities are computed through
% Bayesian theorem, the computed values do not change.

%% 3.b.1 N=10000 bits (Manual Computation)
N_10k = 10000; % number of bits to be generated
tx_10k = rand(1,N_10k) > 0.5; % generate 10000 random bits
tx0_10k = (sum(tx_10k(:)==0))/N_10k % rate of tx 0's from 10k sample bits
tx1_10k = (sum(tx_10k(:)==1))/N_10k % rate of tx 1's from 10k sample bits
txe_10k = 0.05; % transmission error probability rate
txs_10k = 0.95; % transmission success probability rate 

% Manually computed P(R) - vector | manual XOR implementation
manual_rx_10k = tx_10k.*txs_10k+tx_10k.*txe_10k;
%% 

% P(R'1') of 10k bits
rx1_10k = tx1_10k*txs_10k+tx0_10k*txe_10k;
fprintf('Given the conditions above, the probability that the receiver receive 1s from 10000 sample bits is %.4f',rx1_10k);
%% 

% P(R'0') of 10k bits
rx0_10k = tx0_10k*txs_10k+tx1_10k*txe_10k;
fprintf('Given the conditions above, the probability that the receiver receive 0s from 10000 sample bits is %.4f',rx0_10k);
%% 

% P(R'1'|T'1') of 10k bits
PR1T1_10k = (txs_10k*rx1_10k)/tx1_10k;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 1 is %.4f', PR1T1_10k);
%% 

% P(R'0'|T'1') of 10k bits
PR0T1_10k = (txe_10k*rx0_10k)/tx1_10k;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitter transmitted bit 1 is %.4f', PR0T1_10k);
%% 

% P(R'1'|T'0') of 10k bits
PR1T0_10k = (txe_10k*rx1_10k)/tx0_10k;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 0 is %.4f', PR1T0_10k);
%% 

% P(R'0'|T'0') of 10k bits
PR0T0_10k = (txs_10k*rx0_10k)/tx0_10k;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitter transmitted bit 0 is %.4f', PR0T0_10k);
%% 

% P(T'1'|R'1') of 10k bits
PT1R1_10k = (txs_10k*tx1_10k)/((txs_10k*tx1_10k)+(txe_10k*tx0_10k));
fprintf('Given the conditions above, the probability that the transmitter transmitted bit 1 when the received data by the receiver is bit 1 is %.4f', PT1R1_10k);
%% 

% P(T'0'|R'1') of 10k bits
PT0R1_10k = (txe_10k*tx0_10k)/((txe_10k*tx0_10k)+(txs_10k*tx1_10k));
fprintf('Given the conditions above, the probability that the transmitter transmitted bit 0 when then receiver received bit 1 is %.4f', PT0R1_10k);
%% 

% P(T'1'|R'0') of 10k bits
PT1R0_10k = (txe_10k*tx1_10k)/((txe_10k*tx1_10k)+(txs_10k*tx0_10k));
fprintf('Given the conditions above, the probability that the transmitter transmitted bit 1 when the received data by the receiver is bit 0 is %.4f', PT1R0_10k);
%% 

% P(T'0'|R'0') of 10k bits
PT0R0_10k = (txs_10k*tx0_10k)/((txs_10k*tx0_10k)+(txe_10k*tx1_10k));
fprintf('Given the conditions above, the probability that the transmitter send bit 0 when the receiver received bit 0 is %.4f', PT0R0_10k);

%% 3.b.2 N=10000 bits (MATLAB Computation)
%% 

% lines 302-303 must be inserted here!!!
figure(); subplot(4,1,3)
stem(tx_10k);
xlabel('Transmitted 10k bits'); title('MATLAB Computation'); 
% determine is successfully transmitted or not
ch_10k = rand(1,N_10k) > txs_10k;
rx_10k = xor(tx_10k, ch_10k);
%plot the results'
subplot(4,1,4), stem(rx_10k);
xlabel('Received 10k bits');
subplot(4,1,1), stem(tx_10k); xlabel('Transmitted 10k bits'); title('Manual Computation');
subplot(4,1,2), stem(manual_rx_10k); xlabel('Received 10k bits'); 
%% 

% P(R'1') and P(R'0') - MATLAB values

% rate of received 0's from 10k sample bits
matlab_rx0_10k = (sum(rx_10k(:)==0))/N_10k

% rate of received 1's from 10k sample bits
matlab_rx1_10k = (sum(rx_10k(:)==1))/N_10k
%% 

% PR1T1 from MATLAB gen values - 10k bits
matlab_PR1T1_10k = (txs_10k*matlab_rx1_10k)/tx1_10k 
%% 

% PR0T1 from MATLAB gen values - 10k bits
matlab_PR0T1_10k = (txe_10k*matlab_rx0_10k)/tx1_10k 
%% 

% PR0T1 from MATLAB gen values - 10k bits
matlab_PR1T0_10k = (txe_10k*matlab_rx1_10k)/tx0_10k 
%% 

% PR0T0 from MATLAB gen values - 10k bits
matlab_PR0T0_10k = (txs_10k*matlab_rx0_10k)/tx0_10k 
%% 

% PT1R1 from MATLAB gen values - 10k bits
matlab_PT1R1_10k = (txs_10k*tx1_10k)/((txs_10k*tx1_10k)+(txe_10k*tx0_10k)) 
%% 

% PT0R1 from MATLAB gen values - 10k bits
matlab_PT0R1_10k = (txe_10k*tx0_10k)/((txe_10k*tx0_10k)+(txs_10k*tx1_10k)) 
%% 

% PT1R0 from MATLAB gen values - 10k bits
matlab_PT1R0_10k = (txe_10k*tx1_10k)/((txe_10k*tx1_10k)+(txs_10k*tx0_10k))
%% 

% PT0R0 from MATLAB gen values - 10k bits
matlab_PT0R0_10k = (txs_10k*tx0_10k)/((txs_10k*tx0_10k)+(txe_10k*tx1_10k))

%% Discussion for 3b1-2
% * The calculated P(R'1') and P(R'0') values are defined by the
% variables 'rx1_10k' and 'rx0_10k' for manual computation (see sec.3b1)
% and 'matlab_rx0_10k' and 'matlab_rx1_10k'. Although the small variance
% between the manually computed and MATLAB-generated received bits are not visible from the figure above, the
% numbers of these variables show that the there exist an error rate of up
% to ~5% (depending on the computed values of N) within the 
% transmission channel.
%
% * Similarly, the manually computed (sec.b1-c4) and MATLAB-generated (lines 401-439) 
% PR0T0, PR0T1, PR1T0, PR1T1, T1R1, R1T1, R1T0, and R0T0 show a tiny variance to none from each other assuring
% that there is only a chance of up to ~5% that the received or transmitted bits are
% erroneous. Moreover, the success transmission/receving rate of correct 1
% or 0 bits ccording to these probabilities is up to ~95%, which is aligned
% to the provided specification above.
%
% * By inspection, the values from the manually computed probabilities of the samples by
% conditional differ when the same means applied onto the MATLAB generated
% bits. In contrast to when these probabilities are computed through
% Bayesian theorem, the computed values do not change.

%% Conclusion for sec.3
% * Regardless of the length of N or bits transmitted in the BSC
% configuration above, both the manual computations and MATLAB-generated
% simulations and probabilities suggested that the success rate of
% transmission and receiving over the BSC is ~95% with an error rate of
% ~5%. By inspection, it is also observed that the computed probabilities
% by 'conditional method' provide values with tiny variances between the
% manual computations and MATLAB-generated whilst there are none when using
% 'Bayesian method'. Overall, like the corss-validation formulas of confusion matrix, the equations above determine only the
% probabilities of success and error rates of data transmission and receiving the generated
% N-bits over the BSC and the combination of both values of Tx and Rx.
% Moreover, the computed values of these probabilistic equations proved
% that the success Tx(1,0,1,0)/Rx(1,1,0,0) rates of the considered BSC are
% ~95% when alike bits and ~5% when different bit value.

%% II. Non-Symmetric Channel
% Assume the channel is non-symmetric. The probability that a '1' is
% correctly received as '1' is 0.95 while the probability that a '0' is
% correctly received as '0' is 0.85.

%%
% Unlike the previous cases, the transmission rates of 0 and 1 are required
% to be identified since sucessful receive rate of 0 and 1 are already
% provided. Moreover, P(T'1') and P(T'0') must be determined to compute the 
% same probabilities above. The following equations below are necessary to 
% suffice these requirements.

%% 1.a Manual Computation
% * For this section of the lab activity, the author considered generating
% 100 random bits to 

N_100_NS = 100; % number of bits to be generated
rx_100_NS = rand(1,N_100_NS) > 0.5; % generate 100 random bits
rx0_100_NS = (sum(rx_100_NS(:)==0))/N_100_NS % rate of rx 0's from 100 sample bits
rx1_100_NS = (sum(rx_100_NS(:)==1))/N_100_NS % rate of rx 1's from 100 sample bits
rxe1_100_NS = 0.05; % rx 1 error probability rate
rxs1_100_NS = 0.95; % rx 1 success probability rate
rxe0_100_NS = 0.15; % rx 0 error probability rated
rxs0_100_NS = 0.85; % rx 0 success probability rate

%% a.1 P(T'1'): P() txing 1s
% 
% $$P(T'1')=P(R'1')\times P(T'1'\mid R'1')+P(R'0')\times P(T'1'\mid R'0')$$
%

manual_tx1_100_NS = rx_100_NS.*rxs1_100_NS+rx_100_NS.*rxe1_100_NS;
manual_tx0_100_NS = rx_100_NS.*rxs0_100_NS+rx_100_NS.*rxe0_100_NS;
manual_tx_100_NS = manual_tx1_100_NS .* manual_tx0_100_NS; %anding ops

% computes the P(T'1')
tx1_100_NS = rx1_100_NS*rxs1_100_NS+rx0_100_NS*rxe1_100_NS; 
fprintf('Given the conditions above, the probability that the transmitter transmitted 1s from 100 sample bits is %.4f',tx1_100_NS); % display P(T'1')

%% a.2 P(T'0'): P() txing 1s
% 
% $$P(T'0')=P(R'0')\times P(T'0'\mid R'0')+P(R'1')\times P(T'0'\mid R'1')$$
%

% computes the P(T'0')
tx0_100_NS = rxs1_100_NS*rxs0_100_NS+rx1_100_NS*rxe0_100_NS; 
fprintf('Given the conditions above, the probability that the transmitter transmitted 0s from 100 sample bits is %.4f',tx0_100_NS); % display P(T'0')

%% b. PR1T1, PR0T1, PR1T0, PR0T0

%% 

% PR1T1
PR1T1_100_NS = (rxs1_100_NS*rx1_100_NS)/tx1_100_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 1 is %.4f', PR1T1_100_NS);
%% 

% PR0T1
PR0T1_100_NS = (rxe1_100_NS*rx0_100_NS)/tx1_100_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitted data by the transmitter is bit 1 is %.4f', PR0T1_100_NS);
%% 

% PR1T0
PR1T0_100_NS = (rxe0_100_NS*rx1_100_NS)/tx0_100_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 0 is %.4f', PR1T0_100_NS);
%% 

% PR0T0
PR0T0_100_NS = (rxs0_100_NS*rx0_100_NS)/tx0_100_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitted data by the transmitter is bit 0 is %.4f', PR0T0_100_NS);
%% 

% PT1R1
PT1R1_100_NS = (rxs1_100_NS*tx1_100_NS)/((rxs1_100_NS*tx1_100_NS)+(rxe0_100_NS*tx0_100_NS)) 
%% 

% PT0R1
PT0R1_100_NS = (rxe0_100_NS*tx0_100_NS)/((rxe0_100_NS*tx0_100_NS)+(rxs1_100_NS*tx1_100_NS)) 
%% 

% PT1R0
PT1R0_100_NS = (rxe1_100_NS*tx1_100_NS)/((rxe1_100_NS*tx1_100_NS)+(rxs0_100_NS*tx0_100_NS))
%% 

% PT0R0
PT0R0_100_NS = (rxs0_100_NS*tx0_100_NS)/((rxs0_100_NS*tx0_100_NS)+(rxe1_100_NS*tx1_100_NS))

%% 1.b MATLAB Computation

figure(); subplot(2,2,2)
stem(rx_100_NS);
xlabel('Received bits'); title('MATLAB Computation'); 
% determine is successfully received or not
ch1_100_NS = rand(1,N_100_NS) > 0.95;
ch0_100_NS = rand(1,N_100_NS) > 0.85;
ch_100_NS = and(ch1_100_NS, ch0_100_NS);
tx_100_NS = xor(rx_100_NS, ch_100_NS);
%plot the results'
subplot(2,2,4), stem(rx_100_NS);
xlabel('Transmitted bits');
subplot(2,2,1), stem(manual_tx_100_NS); xlabel('Received bits'); title('Manual Computation');
subplot(2,2,3), stem(rx_100_NS); xlabel('Transmitted bits'); 

matlab_tx0_100_NS = (sum(tx_100_NS(:)==0))/N_100_NS % rate of tx 0's from 100 sample bits
matlab_tx1_100_NS = (sum(tx_100_NS(:)==1))/N_100_NS % rate of tx 1's from 100 sample bits

%% 

% PR1T1
matlab_PR1T1_100_NS = (rxs1_100_NS*rx1_100_NS)/matlab_tx1_100_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 1 is %.4f', PR1T1_100_NS);
%% 

% PR0T1
matlab_PR0T1_100_NS = (rxe1_100_NS*rx0_100_NS)/matlab_tx1_100_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitted data by the transmitter is bit 1 is %.4f', PR0T1_100_NS);
%% 

% PR1T0
matlab_PR1T0_100_NS = (rxe0_100_NS*rx1_100_NS)/matlab_tx0_100_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 0 is %.4f', PR1T0_100_NS);
%% 

% PR0T0
matlab_PR0T0_100_NS = (rxs0_100_NS*rx0_100_NS)/matlab_tx0_100_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitted data by the transmitter is bit 0 is %.4f', PR0T0_100_NS);
%% 

% PT1R1
matlab_PT1R1_100_NS = (rxs1_100_NS*matlab_tx1_100_NS)/((rxs1_100_NS*matlab_tx1_100_NS)+(rxe0_100_NS*matlab_tx0_100_NS)) 
%% 

% PT0R1
matlab_PT0R1_100_NS = (rxe0_100_NS*matlab_tx0_100_NS)/((rxe0_100_NS*matlab_tx0_100_NS)+(rxs1_100_NS*matlab_tx1_100_NS)) 
%% 

% PT1R0
matlab_PT1R0_100_NS = (rxe1_100_NS*matlab_tx1_100_NS)/((rxe1_100_NS*matlab_tx1_100_NS)+(rxs0_100_NS*matlab_tx0_100_NS))
%% 

% PT0R0
matlab_PT0R0_100_NS = (rxs0_100_NS*matlab_tx0_100_NS)/((rxs0_100_NS*matlab_tx0_100_NS)+(rxe1_100_NS*matlab_tx1_100_NS))

%% Discussion of section II-1a-b
% * Given the specifications above, the probabilities computed manually are
% defined by the variables 'tx0_100_NS', 'tx1_100_NS' for P(T'0') and
% P(T'1'), 'PR1T1_100_NS', 'PR0T1_100_NS', 'PR1T0_100_NS', 'PR0T0_100_NS'
% for P(R'1'|T'1'), P(R'0'|T'1'), P(R'1'|T'0'), and P(R'0'|T'0'),
% respectively. The variables 'PT1R1_100_NS', 'PT1R0_100_NS', 'PT0R1_100_NS',
% 'PT0R0_100_NS' define the P(T'1'|R'1'), P(T'0'|R'1'), P(T'1'|R'0'), and
% P(T'0'|R'0') respectively. 
%
% The probabilities computed through MATLAB are
% defined by the variables 'matlab_tx0_100_NS', 'matlab_tx1_100_NS' for P(T'0') and
% P(T'1'), 'matlab_PR1T1_100_NS', 'matlab_PR0T1_100_NS', 'matlab_PR1T0_100_NS', 'matlab_PR0T0_100_NS'
% for P(R'1'|T'1'), P(R'0'|T'1'), P(R'1'|T'0'), and P(R'0'|T'0'),
% respectively. The variables 'matlab_PT1R1_100_NS', 'matlab_PT1R0_100_NS', 'matlab_PT0R1_100_NS',
% 'matlab_PT0R0_100_NS' define the P(T'1'|R'1'), P(T'0'|R'1'), P(T'1'|R'0'), and
% P(T'0'|R'0') respectively. 
%
% By inspection, the the manually computed and those through MATLAB
% conditional probabilities ahow similar output only that these numbers are
% farther from the requirement of the Asym BC. In contrast with the
% probabilities computed using BAyesian theorem, the values appear to be
% balanced in both Rx and Tx considering the conbitions of 1 and 0 in both
% ends.Furthermore, we may relate the computed probabilities through
% conditional method to the outcom of Binary Erasure Channel where a
% portion of the output had gone missing or erased due to asymmetrical
% design proabilities of Rx and Tx.

%% 2. MATLAB Asym. N=10000
%
N_10k_NS = 10000; % number of bits to be generated
rx_10k_NS = rand(1,N_10k_NS) > 0.5; % generate 100 random bits
rx0_10k_NS = (sum(rx_10k_NS(:)==0))/N_10k_NS % rate of rx 0's from 100 sample bits
rx1_10k_NS = (sum(rx_10k_NS(:)==1))/N_10k_NS % rate of rx 1's from 100 sample bits
rxe1_10k_NS = 0.05; % rx 1 error probability rate
rxs1_10k_NS = 0.95; % rx 1 success probability rate
rxe0_10k_NS = 0.15; % rx 0 error probability rated
rxs0_10k_NS = 0.85; % rx 0 success probability rate

%% a.1 P(T'1'): P() txing 1s
% 
% $$P(T'1')=P(R'1')\times P(T'1'\mid R'1')+P(R'0')\times P(T'1'\mid R'0')$$
%

manual_tx1_10k_NS = rx_10k_NS.*rxs1_10k_NS+rx_10k_NS.*rxe1_10k_NS;
manual_tx0_10k_NS = rx_10k_NS.*rxs0_10k_NS+rx_10k_NS.*rxe0_10k_NS;
manual_tx_10k_NS = manual_tx1_10k_NS .* manual_tx0_10k_NS; %anding ops

% computes the P(T'1')
tx1_10k_NS = rx1_10k_NS*rxs1_10k_NS+rx0_10k_NS*rxe1_10k_NS; 
fprintf('Given the conditions above, the probability that the transmitter transmitted 1s from 100 sample bits is %.4f',tx1_10k_NS); % display P(T'1')

%% a.2 P(T'0'): P() txing 1s
% 
% $$P(T'0')=P(R'0')\times P(T'0'\mid R'0')+P(R'1')\times P(T'0'\mid R'1')$$
%

% computes the P(T'0')
tx0_10k_NS = rxs1_10k_NS*rxs0_10k_NS+rx1_10k_NS*rxe0_10k_NS; 
fprintf('Given the conditions above, the probability that the transmitter transmitted 0s from 100 sample bits is %.4f',tx0_10k_NS); % display P(T'0')

%% b. PR1T1, PR0T1, PR1T0, PR0T0

%% 

% PR1T1
PR1T1_10k_NS = (rxs1_10k_NS*rx1_10k_NS)/tx1_10k_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 1 is %.4f', PR1T1_10k_NS);
%% 

% PR0T1
PR0T1_10k_NS = (rxe1_10k_NS*rx0_10k_NS)/tx1_10k_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitted data by the transmitter is bit 1 is %.4f', PR0T1_10k_NS);
%% 

% PR1T0
PR1T0_10k_NS = (rxe0_10k_NS*rx1_10k_NS)/tx0_10k_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 0 is %.4f', PR1T0_10k_NS);
%% 

% PR0T0
PR0T0_10k_NS = (rxs0_10k_NS*rx0_10k_NS)/tx0_10k_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitted data by the transmitter is bit 0 is %.4f', PR0T0_10k_NS);
%% 

% PT1R1
PT1R1_10k_NS = (rxs1_10k_NS*tx1_10k_NS)/((rxs1_10k_NS*tx1_10k_NS)+(rxe0_10k_NS*tx0_10k_NS)) 
%% 

% PT0R1
PT0R1_10k_NS = (rxe0_10k_NS*tx0_10k_NS)/((rxe0_10k_NS*tx0_10k_NS)+(rxs1_10k_NS*tx1_10k_NS)) 
%% 

% PT1R0
PT1R0_10k_NS = (rxe1_10k_NS*tx1_10k_NS)/((rxe1_10k_NS*tx1_10k_NS)+(rxs0_10k_NS*tx0_10k_NS))
%% 

% PT0R0
PT0R0_10k_NS = (rxs0_10k_NS*tx0_10k_NS)/((rxs0_10k_NS*tx0_10k_NS)+(rxe1_10k_NS*tx1_10k_NS))

%% 2.b MATLAB Computation

figure(); subplot(2,2,2)
stem(rx_10k_NS);
xlabel('Received bits'); title('MATLAB Computation'); 
% determine is successfully received or not
ch1_10k_NS = rand(1,N_10k_NS) > 0.95;
ch0_10k_NS = rand(1,N_10k_NS) > 0.85;
ch_10k_NS = and(ch1_10k_NS, ch0_10k_NS);
tx_10k_NS = xor(rx_10k_NS, ch_10k_NS);
%plot the results'
subplot(2,2,4), stem(rx_10k_NS);
xlabel('Transmitted bits');
subplot(2,2,1), stem(manual_tx_10k_NS); xlabel('Received bits'); title('Manual Computation');
subplot(2,2,3), stem(rx_10k_NS); xlabel('Transmitted bits'); 

matlab_tx0_10k_NS = (sum(tx_10k_NS(:)==0))/N_10k_NS % rate of tx 0's from 10k sample bits
matlab_tx1_10k_NS = (sum(tx_10k_NS(:)==1))/N_10k_NS % rate of tx 1's from 10k sample bits

%% 

% PR1T1
matlab_PR1T1_10k_NS = (rxs1_10k_NS*rx1_10k_NS)/matlab_tx1_10k_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 1 is %.4f', PR1T1_10k_NS);
%% 

% PR0T1
matlab_PR0T1_10k_NS = (rxe1_10k_NS*rx0_10k_NS)/matlab_tx1_10k_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitted data by the transmitter is bit 1 is %.4f', PR0T1_10k_NS);
%% 

% PR1T0
matlab_PR1T0_10k_NS = (rxe0_10k_NS*rx1_10k_NS)/matlab_tx0_10k_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 1 when the transmitted data by the transmitter is bit 0 is %.4f', PR1T0_10k_NS);
%% 

% PR0T0
matlab_PR0T0_10k_NS = (rxs0_10k_NS*rx0_10k_NS)/matlab_tx0_10k_NS;
fprintf('Given the conditions above, the probability that the receiver receive bit 0 when the transmitted data by the transmitter is bit 0 is %.4f', PR0T0_10k_NS);
%% 

% PT1R1
matlab_PT1R1_10k_NS = (rxs1_10k_NS*matlab_tx1_10k_NS)/((rxs1_10k_NS*matlab_tx1_10k_NS)+(rxe0_10k_NS*matlab_tx0_10k_NS)) 
%% 

% PT0R1
matlab_PT0R1_10k_NS = (rxe0_10k_NS*matlab_tx0_10k_NS)/((rxe0_10k_NS*matlab_tx0_10k_NS)+(rxs1_10k_NS*matlab_tx1_10k_NS)) 
%% 

% PT1R0
matlab_PT1R0_10k_NS = (rxe1_10k_NS*matlab_tx1_10k_NS)/((rxe1_10k_NS*matlab_tx1_10k_NS)+(rxs0_10k_NS*matlab_tx0_10k_NS))
%% 

% PT0R0
matlab_PT0R0_10k_NS = (rxs0_10k_NS*matlab_tx0_10k_NS)/((rxs0_10k_NS*matlab_tx0_10k_NS)+(rxe1_10k_NS*matlab_tx1_10k_NS))

%% Discussion of section II-2a-b
% * Given the specifications above, the probabilities computed manually are
% defined by the variables 'tx0_10k_NS', 'tx1_10k_NS' for P(T'0') and
% P(T'1'), 'PR1T1_10k_NS', 'PR0T1_10k_NS', 'PR1T0_10k_NS', 'PR0T0_10k_NS'
% for P(R'1'|T'1'), P(R'0'|T'1'), P(R'1'|T'0'), and P(R'0'|T'0'),
% respectively. The variables 'PT1R1_10k_NS', 'PT1R0_10k_NS', 'PT0R1_10k_NS',
% 'PT0R0_10k_NS' define the P(T'1'|R'1'), P(T'0'|R'1'), P(T'1'|R'0'), and
% P(T'0'|R'0') respectively. 
%
% The probabilities computed through MATLAB are
% defined by the variables 'matlab_tx0_10k_NS', 'matlab_tx1_10k_NS' for P(T'0') and
% P(T'1'), 'matlab_PR1T1_10k_NS', 'matlab_PR0T1_10k_NS', 'matlab_PR1T0_10k_NS', 'matlab_PR0T0_10k_NS'
% for P(R'1'|T'1'), P(R'0'|T'1'), P(R'1'|T'0'), and P(R'0'|T'0'),
% respectively. The variables 'matlab_PT1R1_10k_NS', 'matlab_PT1R0_10k_NS', 'matlab_PT0R1_10k_NS',
% 'matlab_PT0R0_10k_NS' define the P(T'1'|R'1'), P(T'0'|R'1'), P(T'1'|R'0'), and
% P(T'0'|R'0') respectively. 
%
% By inspection, the the manually computed and those through MATLAB
% conditional probabilities ahow similar output only that these numbers are
% farther from the requirement of the Asym BC. In contrast with the
% probabilities computed using BAyesian theorem, the values appear to be
% balanced in both Rx and Tx considering the conbitions of 1 and 0 in both
% ends.Furthermore, we may relate the computed probabilities through
% conditional method to the outcom of Binary Erasure Channel where a
% portion of the output had gone missing or erased due to asymmetrical
% design proabilities of Rx and Tx.

%% 3. When 60% of the tx bits are '1s'.
% * Regardless of the statistics of 1s and 0s of the transmitted bits by
%the transmitter, whether it is biased by 60% for '1s' or not, the probabilities computed
%using conditional or bayesian means will alter as these equations are
%naturally dependent on the P(T'1') and P(T'0') and of the receiver as well
%P(R'1') and P(R'0'). Also, whether the gicen channel is symmetric or not,
%it also expected that the value of the probabilities might change.
%Consider the codes and executions above or try re-executing the sourcefile
%of this lab report and observe the changes it provide for every sampled
%bits by rand(). By these inspections, we may conclude that the
%probabilities studied above are dependent on the probabilities of the bits
%tx and rx by the channel hence, alteration their value might happen when
%these rx and tx probabilities are changed.