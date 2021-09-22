%% Marwin B. Alejo   2020-20221   EE214_Module1-LabEx1
% * Date Performed (d/m/y): 19/09/2021
% * Date Modified (d/m/y): 22/09/2021

%% I. Analog to Digital Conversion
% * a. Analog time and signal
% define "analog" time and signal
t = 0:0.001:1;                  % time in sec
freq = 5;                       % frequency is 5Hz
analog = 2*sin(2*pi*freq*t);    % analog signal

% plot the analog signal
figure; plot(t, analog)         % plot analog signal
xlabel('time, t(sec)')          % x-axis label

%%
% |Discussion: The defined analog signal above contain five sinusoidal cycles over the
% 1-sec timeframe. The above-generated plot of the defined analog signal
% is consistent with the definition of the analog signal. The 'freq'
% variable of the signal dictates the number of cycles the analog signal
% should have whilst the contant value multiplied with the sine function
% determines the amplitude of the signal.|

%%
% * b. ADC: Sampling and Quantization
% sampling
Fs = 25;                        % sampling frequency is 25Hz
n = 0:1/Fs:1;                   % sampling intervals
sampled = 2*sin(2*pi*freq*n)

% plot the samples of the same figure
hold;
stem(n, sampled);                

% quantization
% define 5 quatization levels as positive and negative integers
quantized = round(sampled)
stairs(n, quantized)

%%
% |Discussion: From the generated plot of the digitized signal above, there are
% five samples per sinusoidal. This is expected since the sampling rate is 1/25.
% The sampling rate when multiplied by 5 is equal to 0.2 sec timeframe which is 
% equal to 1-sec timeframe when multiplied to the number of cycles.|

%%
% * c. ADC Percentage Error

% Computes the percentage error of the actual over quantized values.
PercentageError = (abs(sampled-quantized)./abs(sampled))*100;
PercentageError(1,2:26) % shows the percentage error

%%
% |Discussion: The percentage error of the computed sampled and quantized values
% yields a mean of 7.8518% while individual Percentage Error is shown above.|

%%
% * Exercise 1: Hand Calculation (see the attached data.csv for the computed values)
%%
% |This section consider the sampling and quantizing computations over 1-sec
% timeframe. The actual and quantized values in data.csv are already in the
% absolute form. To determine the sampled(actual) and quantized values, data.csv
% must be loaded first. The author of this report uses MSExcel to manually compute
% each values and verify its similarity with the values computed through MATLAB.|

data = importfile("data.csv", [2, 26]);
%%
% * Sampled(Actual) Value/sample(n) are shown below:
ex1_sampled = data(1:25,1); % load the sampled(actual) values
ex1_sampled = table2array(ex1_sampled); % convert table to vector
ex1_sampled = transpose(ex1_sampled); % transpose the vector
ex1_sampled % display the sampled values

%%
% * Quantized Value/sample(n) are shown below:
ex1_quantized = data(1:25,2); % load the sampled(actual) values
ex1_quantized = table2array(ex1_quantized); % convert table to vector
ex1_quantized = transpose(ex1_quantized); % transpose the vector
ex1_quantized % display the quantized values

%%
% * Corresponding Percentage Error
ex1_PercentageError = (abs(ex1_sampled-ex1_quantized)./abs(ex1_sampled))*100;
ex1_PercentageError % shows the percentage error for exercise 1
%%
% |Discussion: The corresponsing percentage error yields a mean of 7.8530% 
% while individual Percentage Error is shown above.|

%%
% * Exercise 2-3(bonus): Plot the percentage error of the ADC.

figure; plot(n(1,2:26),PercentageError(1,2:26)); title('%Err of ADC in section c.');
figure; plot(n(1,2:26),ex1_PercentageError); title('%Err of ADC in Exercise 1.');

%%
% |Discussion: The ADC percentage error in section c and exercise 1 shows similar results
% except for every 5th sample of the cycles which peaks to 100 in section c
% and cut-off in exercise 1. In general, the two plots show similar results 
% on samples 1 to 4 of every cycle and differ on every 5th sample. The 5th 
% sample of every cycle contributes to the percentage error of the digitized
% analog signal.|

%% II. Matrix Operations
% * 1. Multiplication and Element-wise Multiplication
M1_M = [1 2 3; 4 5 6; 7 8 9]; M1_N = [1 1 1; 1 -8 1; 1 1 1];
M1_A = M1_M*M1_N; % computes the product of M1_M and M1_N matrices
M1_B = M1_M.*M1_N; % computes the element-wise product of M1_M and M1_N
M1_A % shows the product of M1_M and M1_N
M1_B % shows the element-wise product of M1_M and M1_N

%%
% |'*' or 'matrix multiplication' output the linear product of the matrices.
% The number of the rows and columns of the matrices multiplied in this operator must be equal.
% '.*' or 'Element-wise Multiplication' output the element-by-element product
% of the multiplied matrices.| 

%%
% * 2. Solving the unknowns of the Linear Systems

% input values of the linear systems
M2_A = [ 1 -2 3; -1 3 -1; 2 -5 5]; 

% output values of the linear systems
M2_B = [9; -6; 17];

% solve the unknowns of the linear systems
M2_xyz = linsolve(M2_A,M2_B); 
M2_xyz % shows the unknowns of the linear systems

%%
% |The solutions to the unknowns of the linear systems are shown by the variable 'M2_xyz', x=1, y=-1, z=2.| 