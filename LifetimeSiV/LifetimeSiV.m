%
%,------------------------------------------------------------------------,
%| Lifetime fitting for SiV centers        (Lukas Hunold @30/01/21)  V3.7 |
%'------------------------------------------------------------------------'
%
%This program can be used to plot and fit a lifetime measurement. You have
%to provide a dat file including a list with the histogram counts. 
%The program can also do a background correction if you provide a similar
%background file and a deconvolution with the IRF, if you provide a similar
%file for the IRF measurement. A multiexponential fit is also possible.

close all
clear variables


%% ---- Optional functions ------------------------------------------------

%First you have to provide the name of your signal file:
signalFile = 'test';
%and the width of your time bins in ns(!):
BIN_TIME = 16e-3;

%If you have a background file and want to correct with that, activate:
performBgCorrection = false;
%And give here the name of your background file:
bgFile = 'in test not available';
%In this case you can also give a scaling factor for the bg:
BG_SCALING = 1;
%You can also perform a fit of the background to understand it better:
performBgFit = false;

%If you want to plot your lifetime histogram, activate:
showHistogram = false;
%And here you can set the limits for the histogram:
PLOT_START_TIME = 0;
PLOT_END_TIME   = 45;

%If you want to perform a lifetime fit, activate:
performFit = true;
%And give the starting and stop time for the fit here:
FIT_START_TIME  = 0;
FIT_END_TIME   = 5;
%You can have a multiple exponential fit (up to 3), by choosing:
EXP_FIT_ORDER = 1;
%If you use a multidimensional fit, you must provide starting params and 
%bounds later in the script, since its not possible to guess them a priori.
%For a single exp fit you can provide the guessed lifetime here (in ns):
TAU_GUESS = 0.7;

%If you want to deconvolute instrument response function (IRF), activate:
deconvoluteIrf = true;
%And give here your IRF measurement histogramm (counts/bin):
irfFile = 'test_IRF';
%If you want to see a plot of the IRF together with the signal, activate:
showIRF = true;

%If you want to show the fit in a plot, activate:
plotFit = true;

%If you want to show the residuals of the plot, activate:
plotResiduals = true;
%And give the number of bins for the residual histogram:
n_BIN_RESIDUALS = 25;


%% ---- Data Readout ------------------------------------------------------

cd 'Data'

signal = importdata([signalFile,'.dat']);

if performBgCorrection
    bg = importdata([bgFile,'.dat']); %#ok<*UNRCH>
else
    bg = signal*0;
end

if deconvoluteIrf
    irf = importdata([irfFile,'.dat']);
end

cd ..

%Create a timescale from the bin numbers:
time = (0:length(signal)-1)' * BIN_TIME;


%% ---- Plot Data ---------------------------------------------------------

figure('Position',[642,572,560,420])

subplot(2,2,1)
plot(time-PLOT_START_TIME,signal)
hold on
plot(time-PLOT_START_TIME,bg)
title  ('Data of emitter and bg')
xlabel ('t / ns')
ylabel ('counts / bin')
axis   ([0 PLOT_END_TIME-PLOT_START_TIME 0 max(signal)*1.05])
legend ('Emitter Data','Background Data')
grid on

subplot(2,2,2)
plot(time-PLOT_START_TIME,signal)
hold on
plot(time-PLOT_START_TIME,signal-bg*BG_SCALING)
title  ('bg corrected Data')
xlabel ('t / ns')
ylabel ('counts / bin')
axis   ([0 PLOT_END_TIME-PLOT_START_TIME 0 max(signal)*1.05])
legend ('Emitter Data','bg corrected Data')
grid on

subplot(2,2,3)
semilogy(time-PLOT_START_TIME,signal)
hold on
plot(time-PLOT_START_TIME,bg)
title  ('Data of emitter and bg (log)')
xlabel ('t / ns')
ylabel ('counts / bin')
axis   ([0 PLOT_END_TIME-PLOT_START_TIME 0 max(signal)*1.5])
legend ('Emitter Data','Background Data','location','south')
grid on

subplot(2,2,4)
semilogy(time-PLOT_START_TIME,signal)
hold on
plot(time-PLOT_START_TIME,signal-bg*BG_SCALING)
title  ('bg corrected Data (log)')
xlabel ('t / ns')
ylabel ('counts / bin')
axis   ([0 PLOT_END_TIME-PLOT_START_TIME 0 max(signal)*1.5])
legend ('Emitter Data','bg corrected Data','location','south')
grid on


%% ---- Data Ordering with deconvolution ----------------------------------
    
if deconvoluteIrf
    [maxIrf,iMaxIrf] = max(irf);
    irf = irf(iMaxIrf:end);
    [maxSignal,iMaxSignal] = max(signal);
    signal = signal(iMaxSignal:end);
    if performBgCorrection
        bg = bg(iMaxSignal:length(signal)+iMaxSignal);
        scaledBg = bg * BG_SCALING;
        signal = signal - scaledBg;
    end
    if length(signal)>length(irf)
        signal = signal(1:length(irf));
    else
        irf = irf(1:length(signal));
    end
    signal = signal./max(signal);
    irf = irf./max(irf);
    
    time = time(1:length(signal));
end


%% ---- Fit Background ----------------------------------------------------

if performBgFit && performBgCorrection
    
    bgFitFunc = @(bgAmplitude,bgTau,bgOffset,x) ...
        bgAmplitude.*exp(-x./bgTau)+ bgOffset; 

    [bgFit] = fit(time,bg,bgFitFunc,'StartPoint',[max(bg),max(time),0]);

    bgParams = confint(bgFit);
    bgTau    = bgFit.tau_bg;
    errBgTau = bgParams(2,2)-bgTau;

    figure('Position',[642,59,560,420])
    plot(bgFit,time,bg,'kx')
    title  ('Lifetime (black) and Fit (red)')
    xlabel ('t / ns')
    ylabel ('counts / bin')
    grid on 
    legend off
    str= ['\tau = ',num2str(bgTau,3),'\pm',num2str(errBgTau,1), ' ns'];
    annot= annotation('textbox',[0.5 0.76 0.07 0.02], 'String',{ str},...
            'FitBoxToText','on');
    
end


%% ---- Perform deconvolution ---------------------------------------------

if showIRF 
    figure('Position',[642,59,560,420])
    plot(time,signal,'k-')
    hold on
    plot(time,irf,'r-')
    title  ('Signal (black) and IRF (red)')
    xlabel ('t / ns')
    ylabel ('counts / bin')
    grid on 
    legend off
    hold on
end

if deconvoluteIrf
    
    %For using the Wiener-Filter, you have to determine the SNR.
    %Define here, how many elements in the tail you want to use as noise:
    nTAIL = 200;
    %The tail is then defined as the last nTail elements of the signal:
    tail = signal((length(signal)-nTAIL):length(signal));
    %The SNR will then be calculated by this function:
    SNR = rms(signal)/rms(tail)-1;
    %Now you can subtract the noise to have the real signal:
    signal=signal-mean(tail);
    %And perform the FTs to find the deconvoluted signal afterwards:
    FFT_COEFFICIENT = 2^nextpow2(length(signal));
    signalFT = fft(signal,FFT_COEFFICIENT); 
    irfFT = fft(irf,FFT_COEFFICIENT);
    %Now the actual deconvolution is done by the filterfunction
    signalFT = (signalFT./irfFT).*...
        ((irfFT.*conj(irfFT))./((irfFT.*conj(irfFT))+1/SNR));
    %And the back FT
    signal = (ifft(signalFT)+abs(min(ifft(signalFT))));
    %Normalize signal and cut first 2 entries (have usually a wrong value):
    signal = signal(3:length(time)+2)/max(signal(3:length(time)+2));
end


%% ---- Fit signal --------------------------------------------------------

%Cut the signal (and bg) for fitting:
signal = signal(time<FIT_END_TIME);
signal(time<FIT_START_TIME) = [];
time = time(1:length(signal));

if performBgCorrection
    bg = bg(time<FIT_END_TIME);
    bg(time<FIT_START_TIME) = [];
    scaledBg = bg * BG_SCALING;
    signal = signal - scaledBg;
end

% Use exponential fit with order defined in the beginning:

if performFit

    if EXP_FIT_ORDER == 1
        f = @(ampl1,tau1,offset,x) ...
            ampl1.*exp(-x./tau1)+ offset;   

        [expFit,gof,infoout] = fit(time,signal,f,...
                            'StartPoint',[1,TAU_GUESS,0]);

        expParams = confint(expFit);
        tau1 = expFit.tau1;
        errTau1 = expParams(2,2)-tau1;


    elseif EXP_FIT_ORDER == 2
        f = @(amp11,tau1,amp12,tau2,offset,x) ...
            amp11.*exp(-x./tau1) + ...
            amp12.*exp(-x./tau2) + offset;     

        [expFit,gof,infoout] = fit(time,signal,f,...
                                'StartPoint',[21500,1,5500,3,0], ...
                                'lower',[0,0.4,0,1.3,0], ...
                                'upper',[50000,1.3,50000,50,1000]);

        expParams = confint(Fit);
        tau1 = expFit.tau1;
        tau2 = expFit.tau2;
        errTau1 = expParams(2,2)-tau1;
        errTau2 = expParams(2,4)-tau2;


    elseif  EXP_FIT_ORDER == 3
        f = @(amp11,tau1,amp12,tau2,amp13,tau3,offset,x) ...
            amp11.*exp(-x./tau1) + ... 
            amp12.*exp(-x./tau2) + ...     
            amp13.*exp(-x./tau3) + offset;     

        [expFit,gof,infoout] = fit(time,signal,f,...
                       'StartPoint',[10000, 18.00,0,0.25,0,1.0,0], ...
                       'lower',[0,3.00,0,0.10,0,0.5,0], ...
                       'upper',[50000,25.00,50000,0.80,50000,1.0,1000]);

        expParams = confint(Fit);
        tau1 = expFit.tau1;
        tau2 = expFit.tau2;
        tau3 = expFit.tau3;
        errTau1 = expParams(2,2)-tau1;
        errTau2 = expParams(2,4)-tau2;
        errTau3 = expParams(2,6)-tau3;

    end

end


%% ---- Plot Fit ----------------------------------------------------------

if performFit && plotFit
    figure('Position',[1204,574,560,420])
    plot(expFit,time,signal,'k-')
    title  ('Lifetime (black) and Fit (red)')
    xlabel ('t / ns')
    ylabel ('counts / bin')
    grid on 
    legend off
    str= ['\tau = ', num2str(tau1,3), '\pm', num2str(errTau1,1), ' ns'];
    hAnnot= annotation('textbox',[0.5 0.76 0.07 0.02], 'String',{ str},...
            'FitBoxToText','on');
end


%% ---- Plot Residuals ----------------------------------------------------

if performFit && plotResiduals
    figure('Position',[1209,59,560,420])
    if EXP_FIT_ORDER == 1
    residuals = expFit.ampl1.*exp(-time./expFit.tau1)+expFit.offset-signal;
    elseif EXP_FIT_ORDER == 2
    residuals = expFit.ampl1.*exp(-time./Fit.tau1)+ ... 
                expFit.amp12.*exp(-time./Fit.tau2)+expFit.offset-signal; 
    elseif EXP_FIT_ORDER == 3
    residuals = expFit.amp11.*exp(-time./expFit.tau1)+ ... 
                expFit.amp12.*exp(-time./expFit.tau2)+ ...     
                expFit.amp13.*exp(-time./expFit.tau3)+expFit.offset-signal; 
    end
    minRes = min(residuals);
    maxRes = max(residuals);  
    binEdges = linspace(minRes,maxRes, n_BIN_RESIDUALS+1);
    h = histogram(residuals, binEdges);
    title ('Residuals of the Fit')
    xlabel('fit - data')
    ylabel('counts per bin')
end
  
    
%% ---- Protocol of updates -----------------------------------------------

% 26.08.20 (V3.6):  - Sections introduced.
%
% 30.01.21 (V3.7):  - Warnings removed.

