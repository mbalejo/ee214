%% Marwin B. Alejo 2020-20221 EE214_Module3-LabEx1
% *

%% I. Simulation of Random Processes and Stationarity
% *1.a. Stationarity* 
numFxn = 10; % number of sample fxn
samples = 1000; % number of samples/fxn
x = randn(samples,numFxn); % generate ten sample fxns w/ u=0 o=1
t = 0:999; % sample interval

timeMean = mean(x,2);
sampleMean = mean(x',1);

timeVariance = var(x,1,2);
sampleVariance = var(x',1,1);

disp(['1st Moment Stationary- ', num2str(mean([timeMean(:) - sampleMean(:)] .^ 2))]);
disp(['2nd Moment Stationary - ', num2str(mean([timeVariance(:) - sampleVariance(:)] .^ 2))]);

% plot the sample fxns
figure; plot(t,x); title('Fig.1: 10 Sample Fxns (Random Process)');

% plot ensemble average and time average
figure; plot(t,sampleMean); title('Fig.2: Ensemble Mean');
figure; plot(t,timeMean); title('Fig.3: Time Mean');

%%
% *Discussion:* 
% Without further processing and simulations and by just looking with fig1-3, 
% the generated random processes above may be classified as
% Stationary. Stationary process is a stochastic process with its
% statistical propoerties don not alter with time. Its joint probability
% distribution or first/seconds moments are constant. In the case of our
% given vectors, both the means and variance are constant (independent of
% time) and as proven by the yielded output '0' of the first and second
% moments. First and second moments indicate that the vectors (processes)
% differ by elements but not their mean and variance in terms of time.

%%
% *1.b. Ergodic* 

numFxn = 1000; % number of sample fxn
samples = 1000; % number of samples/fxn
x = randn(samples,numFxn); % generate ten sample fxns w/ u=0 o=1
t = 0:999; % sample interval

timeMean = mean(x,2);
sampleMean = mean(x,1);

timeVariance = var(x,1,2);
sampleVariance = var(x,1,1);

disp(['1st Moment Ergodicy- ', num2str(mean([timeMean(:) - sampleMean(:)] .^ 2))]);
disp(['2nd Moment Ergodicy - ', num2str(mean([timeVariance(:) - sampleVariance(:)] .^ 2))]);

% plot the sample fxns
figure; plot(t,x); title('Fig.4: 1k Sample Fxns (Random Process)');

% plot ensemble average and time average
figure; plot(t,sampleMean); title('Fig.5: Ensemble Mean');
figure; plot(t,timeMean); title('Fig.6: Time Mean');

%%
% *Discussion:* 
% Similar to 1a, the vectors produced in 1b may be classified as both
% Stationary (by mere plot observation) and Ergodic as the moment along the
% time and sample dimensions are almost equivalent (see first and second
% moments of ergodicy). Unlike in Stationary that the moment are 'the same
% as is', in Ergodic, the moments are almost or equivalently the same but
% not as is, as always. Moreover, a sample is used to represent the
% entirety of the whole processes. Take note that these moments differ
% for every time execution as they are taken from particular samples of the
% processes and not the entire processes.

%%
% *2a-b. $x(t)=2sin(2*\pi*0.002t)$* 
N = randn(1000,1000); 
t = 0:999; % sample interval
x = 2*sin(2*pi*0.002*t)+N;

figure; plot(t,x(:,1:5)); title('Fig.7: First five sample fxns');

timeMean = mean(x,2);
sampleMean = mean(x,1);

%%
% Time mean
figure; plot(t,timeMean); title('Fig.7a: Time Mean')

%%
% Sample mean
figure; plot(t,sampleMean); title('Fig.7b: Sample Mean')

%%
%

timeVariance = var(x,1,2);
sampleVariance = var(x,1,1);

disp(['1st Moment - ', num2str(mean([timeMean(:) - sampleMean(:)] .^ 2))]);
disp(['2nd Moment - ', num2str(mean([timeVariance(:) - sampleVariance(:)] .^ 2))]);

%%
% *Discussion:* 
% By mere observation of the given expression and not Fig.7, it can be
% implied that the processes are neither stationary nor ergodic for x(t)
% rely on time t. It is not stationary for its variance and mean are not
% constant and time-variant. It is not ergodic for its time mean is not
% equal to or tend to the ensemble average. Observe the values generated
% with first and second moments and their difference in values.
% Furthermore, it is observable that the sample mean and variance mean are
% equivalent to each other.

%%
% *3. Stationarity of randomly phased sinusoids.* 
Fs = 8000;
F1 = 100;
t = 0:1/Fs:5/F1-(1/Fs);
A = 1; K = 10;

% Generate K realizations of the phase
Ph = -pi+(2*pi).*rand(1,K);

% Generate freq mtx
Fmat = t'*F1*ones(1,K);

% Generate the phase mtx
Phmat = ones(size(t'))*Ph;

% Generate the sinusoid
Snk = A*sin(2*pi*Fmat+Phmat);

% Plot K relaizations (different phase) of the sinusoids
figure; plot(Snk); title('Fig.8: K realizations of sine');

% Visualize the random process in 3D
figure; mesh(Snk); title('Fig.9: Random Processes');

% compute the time average of each realization
Tave = mean(Snk);

% Ensemble mean
Eave = mean(Snk');

% plot autocorrelation
figure; plot(xcorr(Snk)); title('Fig.10: Autocorrelation of Snk')

timeMean = mean(Snk,2);
sampleMean = mean(Snk',1);

timeVariance = var(Snk,1,2);
sampleVariance = var(Snk',1,1);

disp(['1st Moment WSS- ', num2str(mean([timeMean(:) - sampleMean(:)] .^ 2))]);
disp(['2nd Moment WSS- ', num2str(mean([timeVariance(:) - sampleVariance(:)] .^ 2))]);

%%
% *Proof of non-ergodic:* if we increase K to 100...
Fs = 8000;
F1 = 100;
t = 0:1/Fs:5/F1-(1/Fs);
A = 1; K = 100;
Ph = -pi+(2*pi).*rand(1,K);
Fmat = t'*F1*ones(1,K);
Phmat = ones(size(t'))*Ph;
Snk = A*sin(2*pi*Fmat+Phmat);
Tave = mean(Snk);
Eave = mean(Snk');
timeMean = mean(Snk,2);
sampleMean = mean(Snk',1);
timeVariance = var(Snk,1,2);
sampleVariance = var(Snk',1,1);
disp(['1st Moment Ergodicy Test- ', num2str(mean([timeMean(:) - sampleMean(:)] .^ 2))]);
disp(['2nd Moment Ergodicy Test- ', num2str(mean([timeVariance(:) - sampleVariance(:)] .^ 2))]);

%%
% *a-b.* 
% In a generalized statement, the randomly geenrated sinusoids above may
% be classified as a Wide-Sense Stationary (WSS) random process. The
% reasoning for this claim is similar to the general description of a
% stationary random process in section 1 above. Moreover, it is evident
% that the generated sinusoidal processes are time-invariant with its
% moments yield a constant '0' for the first and second moment test (see
% 1st/2nd moment WSS above). 
%
% The randomly generated sinusoid may also be classified as an
% Ergodic Random Process. Ofcourse, to suffice the definition of
% ergodicity, consider an increased realization K to 100. Generating the
% sinusoids time mean and sample mean yield approximately equivalent results as shown in
% the _proof-section above_. Without further words and from this
% observation, it can be inferred that this random process is an Ergodic RP.
%
% Given Fig.10 or the autocorrelation plot of the generated random process
% with 10 realizations, aside from the generic stationary properties, it
% suffices the ff. WSS properties among others: *(1)* $R_{x}(t)$ is an even 
% function given $R_{x}(t) = R_{x}(-t)$ if we consider the x=400 of fig.10 
% as our time origin. *(2)* It is clear that the generated random process is 
% a WSS given that its autocorrelation plot achieve a maximum realization 
% at time origin (t=0) which is 400 in fig.10.

%%
% *Discussion:* 
% Regardless if the random processes are generated using sin/cos fxns or
% randn/rand, a random process is stationaru if *(1)* its mean is a contant
% (time-invariant) value $u_{x}(t)=ux$; *(2)* its mean square value is also a
% constant value; *(3)* its variance is constant value
% $\sigma_{x}^{2}(t)=\sigma_{x}^{2}$; *(4)* its autocorrelation depends on
% the time distance between two samples $R_{x}(t_{1},t_{2})=R_{x}(T)$; *(5)*
% its autocovariance depends on time distance between two samples at time
%  $t_{1}$ and $t_{2}: K_{x}(t_{1},t_{2})=K_{x}(T)$ -- _be noted that
% autocovariance and autocorrelation depends on time but the entire
% functions of the process itself_ .

%%
% *4. Real-life realization of Random Process* 

%%
% *a. Paying monthly electricity bill.* 
% Paying monthly electric bill requires an individual to perform certain
% processes according to their preferred mode of payment. Considering the
% traditional electric bill over-the-counter payment, an individual
% requires to fall-in-line prior to his/her bill to be received by the
% counter and the volume of the line of people paying their bill
% over-the-counter depend on the time of the day. In general, the
% randomness of this process (electric bill payment) is the number/volume of
% people paying their monthly electric bill over-the-counter varies each
% time(minutes,hours,days,weeks,monthly). This volume of bill payers may be
% averaged for a set of time and may be used to statistically determine the
% possible volume of electric bill customers who are going to pay through
% over-the-counter payment mode for a particular time.

%%
% *b. Listening or watching the daily Covid-19 news cases.* 
% Similar to 4a above, the randomness of this process are the numbers or
% volume of observed/reported Covid-19 cases per day and when averaged
% by a factor of month, semi-annual, or annual. This averaged number may be
% used as an activation value for determining the possible volume of case
% for a particular time signature through the observed parameters.

%%
% *c. Observations of people walking and vehicles passing in a busy
% intersection.* 
% Similar to 4a and 4b, the volume of people and vehicles passing through
% an intersection differ in all moments of time. Given a day, the number or
% volume of people crossing the pedestrians of the intersection and the
% number of cars passing-through the same intersection differ every hour
% for a total number of 24 random processes for a day. A collection of this
% record for a month might be used to create a predicitve model for the
% later months which would allow the possible numbers of people and cars
% passing through the same intersection in each hour of the following
% months and times -- ofcourse in consideration of the observed parameters.

%%
% *Why is each process a random one?* 
% Although it is considered that in real-life, nothing may be done twice or
% exactly the same under the same time but only a particular process may
% occupy a particular moment of time per function. Furthermore, there are
% processes that may be considered to as stationary or repetitive, the
% values or function performed under a certain still vary from one process
% execution to another hence random.