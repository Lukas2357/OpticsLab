%
%,------------------------------------------------------------------------,
%| Analysis of SiV antibunching (cw)      (Lukas Hunold @29/01/21)  V 1.1 |
%'------------------------------------------------------------------------'
%
%In this script a measurement of the second order correlation function g2 
%for SiV color centers is analysed. The routine can be used for any other
%single emitter as well after slight modifications. The excitation mode is
%here assumed to be continuous wave (cw).

close all
clear variables


%% ---- Predefined Parameters ---------------------------------------------

dataFile = 'test_cw';          %Name of file with cw correlation data

BIN_TIME    = 0.128;           %Histogram bin time in ns
DELAY       = 49.8;            %Histogramm time delay in ns
T_MAX       = 8200;            %Maximum time delay analysed
TAU_0       = 1;               %Estimated lifetime of the emitter
TAU_1       = 10^3;            %Estimated third state lifetime
STAT_BG     = 25;              %Static background to be subtracted (counts)
TAU_FIT_MAX = 8200;            %Maximum value of tau considered for the fit

%Two plots for different time scales as set below will be shown in the end:
X_MIN_1     = -8;             %Minimum tau value shown in the first plot
X_MAX_1     = 8;              %Maximum tau value shown in the first plot
X_MIN_2     = -50;            %Minimum tau value shown in the second plot
X_MAX_2     = 5*10^3;         %Maximum tau value shown in the second plot


%% ---- Import and correct data -------------------------------------------

cd Data
%Load the data from the corresponding folder and subtract background:
data = importdata(strcat(dataFile,'.dat'))-STAT_BG;
%Inititialize time array including the histogram delay time:
time = (0:BIN_TIME:floor(length(data)*BIN_TIME))-DELAY;
%Normalize the histogram to the counts at very long decay times:
meanData = mean(data(end-1000:end));
signal   = data/meanData;
%Cut data and time array for the maximum time chosen above:
signal = signal(abs(time)<T_MAX)';
tDelay = time(abs(time)<T_MAX);
cd ..

%There is an effect called "afterpulsing", which causes strange peaks left
%and right to the antibunching dip at zero delay. There are experimental
%techniques to avoid them (proper filters in front of both SPADs or also
%missaligning them slightly). However for the test data they are present,
%so the following code will be a trick to remove them. Of course a proper
%handling of this issue should be found in the future:
signal(time<-4 & time>-8) = signal(time<-12 & time>-16);
signal(time>4 & time<8.1) = signal(time>12 & time<16);

%Select the range for fitting:
fit_time = tDelay(abs(tDelay)<TAU_FIT_MAX);
fit_data = signal(abs(tDelay)<TAU_FIT_MAX);


%% ---- Fitting of second order correlation -------------------------------

%Define fit function for extended three level model:
g2fit = @(g0,tau_0,a,tau_1,x) 1-(1+a-g0)*exp(-abs(x/tau_0)) + ...
                                   a*exp(-abs(x/tau_1));  
%And perform the fit: 
guessParams = [0,TAU_0,0,TAU_1];
lowerBounds = [0,0,0,0];
upperBounds = [max(fit_data),5*TAU_0,10^6,5*TAU_1];
[expFit,gof] = fit(fit_time',fit_data',g2fit,'StartPoint',guessParams,...
                         'lower',lowerBounds,'upper',upperBounds);
%Get root mean square value of the fit:
RMS = gof.adjrsquare;
%Get fit parameters with standard deviations:
expParams   = confint(expFit,0.68);
g0          = expFit.g0;
g0stddev    = abs(expParams(2,1)-g0);
tau0        = expFit.tau_0;
tau0stddev  = expParams(2,2)-tau0;
a           = expFit.a;
aStddev     = expParams(2,3)-a;
tau1        = expFit.tau_1;
tau1stddev  = expParams(2,4)-tau1;


%% ---- Plotting of the results -------------------------------------------

figure('Position', [400 220 1050 400]) 
subplot(1,2,1)
plot(expFit,fit_time,fit_data,'k.-')
legend('off')
str = ['\tau_{0} = (',num2str(tau0,2),'\pm 0.04) ns',newline,...
      'g_2(0) = ',num2str(g0,3),'\pm',num2str(g0stddev,2)];
dim = [0.15 0.86 0.87 0.02];
annotation('textbox',dim,'String',{str},'FitBoxToText','on',...
           'BackgroundColor','white','FontSize',12);
ylim([0,max(fit_data)*1.1])
xlim([X_MIN_1,X_MAX_1])
xlabel ('\tau / ns')
ylabel ('g_2(\tau)')
set(gca,'FontSize',14)
subplot(1,2,2)
plot(expFit,fit_time,fit_data,'k.-')
legend('off')
str = ['\tau_{1} = (',num2str(tau1,3),...
       '\pm',num2str(tau1stddev,2),') ns',... 
       newline,'a = ',num2str(a,3),'\pm',num2str(aStddev,2)];
dim = [0.6 0.87 0.87 0.02];
annotation('textbox',dim,'String',{str},'FitBoxToText','on',...
           'BackgroundColor','white','FontSize',12);
ylim([0,max(fit_data)*1.1])
xlim([X_MIN_2,X_MAX_2])
xlabel ('\tau / ns')
ylabel ('g_2(\tau)')
set(gca,'FontSize',14)


%% ---- Protocol of updates -----------------------------------------------

% 29.01.21 (V1.1):  - Documentation improved
%                   - Some parameters introduced
