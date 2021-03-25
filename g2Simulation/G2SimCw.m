%
%,------------------------------------------------------------------------,
%| Simulation of SiV antibunching (cw)    (Lukas Hunold @01/02/21)  V 1.2 |
%'------------------------------------------------------------------------'
%
%In this script a measurement of the second order correlation function g2 
%for a single emitter in continuous mode is simulated.

close all
clear variables
%#ok<*NASGU>
%#ok<*UNRCH>


%% ---- Optional functions ------------------------------------------------

%If you want to display the plots of the g2 histogram, activate:
showPlots = 1;
%And set the waiting time after any of them here to watch them properly:
PLOT_DELAY = 1;
%If you want to normalize the g2 function, activate:
normG2 = 0;
%If you want to see the progress of the simulation, activate:
printStatus = 1;

%You can vary different parameters for the subsequent simulations in order
%to check their effect on the g2 outcome and the fit results. Here you can
%decide which parameter should be varied (make sure to always activate one 
%option only, otherwise errors will occur!):

%Change measurement time (~total collected signal) in each sim:
varyMeasureTime = 0;
%Change total cps on SPAD in each sim (signal/bg ratio kept constant):
varyTotalCps    = 1;
%Change only signal cps in each sim (signal/bg ratio increases):
varySignalCps   = 0;         

%If you activate one of these, a number of nREP simulations will be
%performed, changing the corresponding parameter in equal steps up to the
%value you set below. If you deacivate it, only the value set below is
%used. Set the number of simulation repetitions nREP here:
nREP     = 10;             %Repetitions of the simulation
%If you choose your parameters such that there are by far too few events, a
%constant default value will be used instead. So if your parameters do not
%vary, although you told them to do so, just increase any of the values.


%% ---- Input Parameters --------------------------------------------------

%These three are the parameters corresponding to the three options above:
T_MEASURE = 36000;              %Measurement time in seconds
TOT_CPS   = 10000;              %Mean recorded counts per second
SIG_CPS   = 10000;              %Signal counts out of this (rest bg)

%These are the other relevant parameters of the simulated experiment:
TAU_1    = 1;                          %Lifetime of the emitter in ns
APD_TIME = 0.004;                      %Fastest recording time of the APDs
BIN_TIME = [0.016,0.064,0.256,1.024];  %Bin time length of histogram in ns
TAU_MAX  = 4;                          %Maximum recorded time delay in ns
G0       = 0;                          %Offset value of G2 at tau=0

%The following value will be used to guarantee, that you have at least one
%entry in the histogram, since the normalization will otherwise diverge. If
%you see no data points in your plot and get an error, increase this value:
MIN_CRIT = 10^6;


%% ---- Setting up --------------------------------------------------------

%First check if by accident the signal cps are higher than the total cps
%and set them equal in that case to avoid unphysical results:
if TOT_CPS < SIG_CPS
    TOT_CPS = SIG_CPS;
end
%Then initialize time array, g2-fct array for the ideal TLS, bin time 
%array, fit parameter arrays and arrays for varying params:
tau_raw  = -TAU_MAX:APD_TIME:TAU_MAX;
g2_raw   = 1-(1-G0)*exp(-abs(tau_raw/TAU_1));
binTime  = zeros(4,1);
g0       = zeros(4,nREP);
g0std    = zeros(4,nREP);
tau1     = zeros(4,nREP);
tau1std  = zeros(4,nREP);
RMS      = zeros(4,nREP);
tMeasure = ones(nREP,1)*T_MEASURE;
totCps   = ones(nREP,1)*TOT_CPS; 
sigCps   = ones(nREP,1)*SIG_CPS;
%The three parameters are only varied if the option is activated. In
%addition it is guaranteed here, that a minimum count per bin is present to
%avoid divergence of the normalization:
if varyMeasureTime
    minTime   = max(10^2*MIN_CRIT/SIG_CPS^2,T_MEASURE/nREP); 
    maxTime   = max(10^2*MIN_CRIT/SIG_CPS^2,T_MEASURE);
	tMeasure  = linspace(minTime,maxTime,nREP);
end
if varyTotalCps
    minTotCps = max(MIN_CRIT/T_MEASURE*(SIG_CPS/TOT_CPS)^2,TOT_CPS/nREP);
    maxTotCps = max(MIN_CRIT/T_MEASURE*(SIG_CPS/TOT_CPS)^2,TOT_CPS);
	totCps    = linspace(minTotCps,maxTotCps,nREP);
    sigCps    = totCps*SIG_CPS/TOT_CPS;
end
if varySignalCps
    minSigCps = max(MIN_CRIT/T_MEASURE,SIG_CPS/nREP);
    maxSigCps = max(MIN_CRIT/T_MEASURE,SIG_CPS);
	sigCps    = linspace(minSigCps,maxSigCps,nREP);
end
%Next calculate background by difference of total and signal counts:
bgCps = totCps-sigCps;


%% ---- Simulating and fitting of g2 measurements -------------------------

%Open figure to display simulated g2 functions and fits:
if showPlots
    figure('Position', [220 120 1400 820]) 
end
%And start the simulation:
for iREP = 1:nREP
    
    %Initialize cell arrays for bunch of g2 functions, time axes and fits:
    g2s    = cell([1 4]);
    taus   = cell([1 4]);
    expFit = cell([1 4]);
    
    %Iterate over different binnings of the histogram:
    for iBinning = 1:4
        %Vary here the binning time according to your settings:
        if length(BIN_TIME) == 4
            binTime(iBinning) = BIN_TIME(iBinning);
        else
            binTime(iBinning) = 0.256;
            inputError = strcat('BIN_TIME must be a 4-entry array ',...
                              ' using default bin time of 256ps instead');
            if (iREP == 1 || iREP == nREP) && iBinning == 4
                fprintf(inputError)
            end
        end
        
        %Calculate according bin parameters and initialize g2/tau arrays:
        bin_length = floor(binTime(iBinning)/APD_TIME);
        n_bins     = floor(length(g2_raw)/bin_length);
        g2         = zeros(1,n_bins);
        tau        = zeros(1,n_bins);
        
        %Perform binning based on above parameters:
        for iTau=1:n_bins
            bin       = g2_raw((iTau-1)*bin_length+1:iTau*bin_length);
            g2(iTau)  = mean(bin);
            tau(iTau) = tau_raw(iTau*bin_length);
        end
        taus{iBinning} = tau;
        
        %Get the average bin counts (for signal and background):
        bin_counts   = totCps(iREP)^2*...
                        binTime(iBinning)*10^-9*tMeasure(iREP);
        s_bin_counts = (totCps(iREP)-bgCps(iREP))^2*...
                        binTime(iBinning)*10^-9*tMeasure(iREP);
        b_bin_counts = bin_counts - s_bin_counts;
        
        %Add Poissonian noise to the g2s and normalize if activated:
        g2s{iBinning}  = g2.*poissrnd(s_bin_counts,size(g2)) + ...
                             poissrnd(b_bin_counts,size(g2));
        if normG2
            g2s{iBinning}  = round(g2s{iBinning})./ceil(bin_counts);
            normFactor = 1;
        else
            normFactor = bin_counts;
        end
        
        %Initialize fit fct, normalized fit data array and bounds:
        fitfct = @(g0,tau1,scale,x) (1-(1-g0)*exp(-abs(x/tau1)))*scale;
        lb     = [0,0,normFactor];
        ub     = [max(g2s{iBinning}),3*TAU_1,normFactor];
        guess  = [G0,TAU_1,normFactor];
        
        %Perform the fit:
        [expFit{iBinning},gof] = fit(tau',g2s{iBinning}',fitfct,...
                                 'StartPoint',guess,'lower',lb,'upper',ub);
                                
        %Get RMS value and fit parameters with standard deviations:                         
        RMS(iBinning,iREP)     = gof.adjrsquare;
        expParams              = confint(expFit{iBinning},0.68);
        g0(iBinning,iREP)      = expFit{iBinning}.g0;
        g0std(iBinning,iREP)   = abs(expParams(2,1)-g0(iBinning,iREP));
        tau1(iBinning,iREP)    = expFit{iBinning}.tau1;
        tau1std(iBinning,iREP) = expParams(2,2)-tau1(iBinning,iREP);
        
    end
    
    %Restart the loop for faster acquiring of subsequent plots:
    if showPlots
        for iBinning = 1:4
            subplot(2,2,iBinning)
            %Plot randomized g2 function and exp. fit:
            plot(expFit{iBinning},taus{iBinning},g2s{iBinning},'k-')
            s1 = sprintf('total cps: %.0f ', totCps(iREP));
            s2 = sprintf(' // signal cps: %.0f ', sigCps(iREP));
            s3 = sprintf(' // integration time: %.0f min %.0f s ',...
                         tMeasure(iREP)/60,mod(tMeasure(iREP),60));
            title(strcat(s1,s2,s3))
            xlabel ('\tau / ns')
            ylabel ('g_2(\tau)')
            xlim([-TAU_MAX,TAU_MAX])
            ylim([0,max(g2s{iBinning})*1.3])
            legend off
            hold on
            %Display tau_1 and g2(0) with errors:
            str1 = ['\tau = ',num2str(tau1(iBinning,iREP),2), '\pm', ...
                              num2str(tau1std(iBinning,iREP),2),' ns', ... 
                    newline, ...
                    'g_2(0) = ',num2str(g0(iBinning,iREP),2), '\pm', ...
                                num2str(g0std(iBinning,iREP),1)];
            dim1 = {[0.15 0.875 0.07 0.02] [0.6 0.875 0.07 0.02]
                   [0.15 0.4 0.07 0.02] [0.6 0.4 0.07 0.02]};
            Annot1 = annotation('textbox',dim1{iBinning},'String',...
                     {str1},'Color','r','FitBoxToText','on',...
                     'BackgroundColor','white', 'EdgeColor','None');
            str2 = sprintf('bin time: %.0f ps',binTime(iBinning)*1000);
            dim2 = {[0.35 0.875 0.07 0.02] [0.8 0.875 0.07 0.02]
                   [0.35 0.4 0.07 0.02] [0.8 0.4 0.07 0.02]};
            Annot2 = annotation('textbox',dim2{iBinning}, 'String',...
                     {str2},'FitBoxToText','on','BackgroundColor','white');
        end
        %Wait for the next iteration to be displayed:
        pause(PLOT_DELAY)
        if iREP ~= nREP
            delete(findall(gcf,'type','line'))
            delete(findall(gcf,'type','annotation'))
        end
    end
    if printStatus
        status = sprintf('repetition %.0f/%.0f done',iREP,nREP);
        fprintf(status)
        fprintf('\n')
    end
end


%% ---- Plot final results ------------------------------------------------

figure('Position', [220 120 1400 820]) 
binLabel = cell(4,1);
for i = 1:4
    binLabel{i} = strcat(num2str(binTime(i)*1000),' ps bin time ');
end
options = [varyTotalCps,varyTotalCps,varyMeasureTime];

plotResult(totCps,sigCps,tMeasure,g0','g_2(0)',binLabel,options,1)
title  ('g_2(0) from fit')

plotResult(totCps,sigCps,tMeasure,g0std','g_2(0) std',binLabel,options,2)
title ('g_2(0) standard deviation')

plotResult(totCps,sigCps,tMeasure,tau1','\tau_1',binLabel,options,3)
title ('\tau_1 from fit')
ylim([0,2*TAU_1])

plotResult(totCps,sigCps,tMeasure,tau1std','\tau_1 std',binLabel,options,4)
title ('\tau_1 standard deviation')
ylim([0,2*TAU_1])

plotResult(totCps,sigCps,tMeasure,RMS','fit RMS',binLabel,options,[5,6])
title ('fit adjusted root mean square')


%% ---- Used Functions ----------------------------------------------------

function plotResult(x1,x2,x3,yAxis,yLabel,binLabel,conditions,nSubplot)
    subplot(3,2,nSubplot)
    if conditions(1)
        plot(x1,yAxis);
        xlabel('mean total cps on SPAD')
    elseif conditions(2)
        plot(x2,yAxis);
        xlabel('mean signal cps on SPAD')
    elseif conditions(3)
        plot(x3/60,yAxis);
        xlabel('measurement time in min')
    else
        plot(yAxis);
        xlabel('repetition of the Simulation')
    end
    ylabel(yLabel);
    ylim([0,1])
    legend(binLabel{1},binLabel{2},binLabel{3},binLabel{4},...
           'Location','Best')
    set(gca,'FontSize',12)
end


%% ---- Protocol of updates -----------------------------------------------

% 30.01.21 (V1.1):  - Introduction of sections
%                   - Proper option for varying the parameters
%                   - Fct for simpler plotting of final results
%                   - Minimum criterion introduced to avoid divergence
%                   - Appearance and documentation improved
%
% 01.02.21 (V1.2):  - Warnings removed
