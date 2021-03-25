%
%,------------------------------------------------------------------------,
%| Analysis of SiV antibunching (pulsed)  (Lukas Hunold @29/01/21)  V 1.1 |
%'------------------------------------------------------------------------'
%
%In this script a measurement of the second order correlation function g2 
%for SiV color centers is analysed. The routine can be used for any other
%single emitter as well after slight modifications. The excitation mode is
%here assumed to be pulsed.

close all
clear variables


%% ---- Predefined Parameters ---------------------------------------------

dataFile = 'test_pulsed';       %Name of file with pulsed correlation data

BIN_TIME     = 0.256;           %Histogram bin time in ns
POST_BINNING = 1;               %Combine this number of bins for analysis
DELAY        = 50;              %Histogramm time delay in ns
T_MAX        = 50;              %Maximum time delay analysed
f0           = 0.08;            %Laser repetition rate in GHz
TAU_0        = 1;               %Estimated lifetime of the background


%% ---- Import and correct data -------------------------------------------

cd Data
%Load the data from the corresponding folder and subtract background:
data_raw = importdata(strcat(dataFile,'.dat'));
time_raw = (0:BIN_TIME:floor(length(data_raw)*BIN_TIME))-DELAY;
data_raw = data_raw(abs(time_raw)<T_MAX)';
time_raw = time_raw(abs(time_raw)<T_MAX);
cd ..

%There is an effect called "afterpulsing", which causes strange peaks left
%and right to the antibunching dip at zero delay. There are experimental
%techniques to avoid them (proper filters in front of both SPADs or also
%missaligning them slightly). However for the test data they are present,
%so the following code will be a trick to remove them. Of course a proper
%handling of this issue should be found in the future:
zero_time_intervall = time_raw(abs(time_raw)<4);
zero_data_intervall = data_raw(abs(time_raw)<4);
data_raw(time_raw<-4 & time_raw>-7.8) = data_raw(time_raw<-4-1/f0 & ...
                                                 time_raw>-8-1/f0);
data_raw(time_raw>4 & time_raw<8) = data_raw(time_raw>4+1/f0 & ...
                                             time_raw<8+1/f0);

%Normalize the data to its maximum:
data_raw = data_raw/max(data_raw);

%Calculate bin parameters and initialize g2/tau arrays:
bin_time = BIN_TIME*POST_BINNING;
n_bins = floor(length(data_raw)/POST_BINNING);
data = zeros(1,n_bins);
time = zeros(1,n_bins);
%Perform binning based on above parameters:
for iTau=1:n_bins
    bin_signal = data_raw((iTau-1)*POST_BINNING+1:iTau*POST_BINNING);
    data(iTau) = mean(bin_signal);
    time(iTau) = time_raw(iTau*POST_BINNING);
end


%% ---- Fitting the measurement result ------------------------------------

%Perform exponential fit of the region around zero delay:
expFct = @(g0,tau_1,x) g0*exp(-abs(x/tau_1));   
time_fit = time(time>-1/f0/2 & time<1/f0/2);
data_fit = data(time>-1/f0/2 & time<1/f0/2);
guessParams = [1,TAU_0];
lowerBounds = [0,0];
upperBounds = [2,5*TAU_0];
[expFit,gof] = fit(time_fit',data_fit',expFct,'StartPoint',guessParams,...
                          'lower',lowerBounds,'upper',upperBounds);
%Get root mean square value of the fit:          
RMS = gof.adjrsquare;
%Get fit parameters with standard deviations:
expParams = confint(expFit,0.68);
g0        = expFit.g0;
%Get the data at zero delay to compare it with the value of the fit:
g0Data    = max([data_raw(time_raw == max(time_raw(time_raw<0))),...
                 data_raw(time_raw == min(time_raw(time_raw>0)))])/...
                 max(data_raw);
g0stddev  = abs(expParams(2,1)-g0);
tau1      = expFit.tau_1;
tau1std   = expParams(2,2)-tau1;


%% ---- Plotting of the results -------------------------------------------

figure('Position', [400 520 650 300])
plot(expFit,time,data,'k-')
legend('off')
hold on
str= ['\tau_{BG} = ',num2str(tau1,2),'\pm',num2str(tau1std,2),' ns', ... 
      newline,'g_2(0) = ',num2str(g0,2),'\pm',num2str(g0stddev,1)];
dim = [0.68 0.875 0.07 0.02];
annotation('textbox',dim,'String',{str},'FitBoxToText','on',...
           'BackgroundColor','white');
ylim([0,max(data)*1.1])
xlabel ('\tau / ns')
ylabel ('g_2(\tau)')
set(gca,'FontSize',11)


%% ---- Protocol of updates -----------------------------------------------

% 29.01.21 (V1.1):  - Documentation improved
%                   - Some parameters introduced
