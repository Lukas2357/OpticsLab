%
%,------------------------------------------------------------------------,
%| Simulation of ensemble saturation      (Lukas Hunold @01/02/21)  V 1.3 |
%'------------------------------------------------------------------------'
%
%In this script the saturation behaviour of an emitter ensemble is
%simulated. We assume to have a number of similar emitters in the confocal
%volume, which are studied at the same time by collecting their emission
%signal for different excitation powers. The saturation powers of the
%individual emitters might be different, which can be taken into account by
%simulating a uniform distribution of P_sat around a mean value with a
%given range. The spatial distribution of emitters is always assumed to be
%uniform and the laser probes a limited area out of it, which results in a
%log-like behaviour of the overall saturation curve, as introduced below.
%
%In addition to that, the effect of different numbers of emitters in the
%confocal volume and of different maximum powers ised to generate the
%curves, can be studied. You will be able to set all the mentioned
%parameters by hand and generate saturation curves together with standard
%saturation fits and log-fits for emitter ensembles. 

close all
clear variables
%#ok<*UNRCH>
%#ok<*NASGU> 


%% ---- Optional functions ------------------------------------------------

%Below you will set a number of emitters and a range of saturation powers
%to be simulated. The saturation curve obtained for the two parameters
%given there will always be generated. If you want to compute also all
%numbers of emitters up to the number of nEmitters set below, activate:
SimEmitterNumbers = 1;
%The same can be done for all ranges of saturation powers up to the one 
%given as P_SAT_RANGE below, by activating
SimPsatRanges     = 0;
%And also for the maximum power MAX_POWER up to which the saturation curve
%is generated, by activating
SimMaxPowers      = 0;
%If you want to show ALL simulated saturation curves with the corresponding
%fits (attention! Based on the inputs below this can be a lot!), activate 
plotSatCurves     = 1;


%% ---- Set initial parameters --------------------------------------------

%First set the characteristic parameters of your emitter ensemble:
P_SAT_MEAN      = 10;       %Mean saturation power of the probed emitters
P_SAT_RANGE     = 20;       %Range of saturation powers around that mean
nEMITTERS       = 20;       %Number of emitters in the confocal volume

%And then the parameters for your saturation measurement:
nPOWER_STEPS    = 10;       %Number of power steps probed for the sat curve
MAX_POWER       = 10;       %Maximum power probed for the sat curves
nRUNS           = 10;        %Number of runs for each given set of params

%And finally the number of steps in which the params are changed in order
%to check their influence on the result. The first one is relevant if
%SimPsatRanges is active, the second if SimMaxPowers is active and the
%third one if SimEmitterNumbers is active:
SAT_RANGE_STEPS = 5;       %Number of steps for ranges up to value above
MAX_POWER_STEPS = 5;       %Number of steps in which max power is variied
nEMITTER_STEPS  = 5;       %Number of steps in which nEmitter varies


%% ---- Initialize arrays -------------------------------------------------

%Initialize power array by nPOWER_STEPS uniform values up to MAX_POWER.
%Increase the number of steps to at least 3 for the lowest maximum power
%simulated, since the fit will otherwise fail due to few data points:
powerSteps = MAX_POWER/max(nPOWER_STEPS,3*MAX_POWER_STEPS);
power      = powerSteps:powerSteps:MAX_POWER;
%Initialize arrays for the emitter ensemble signal at a given power
%("counts") and the saturation powers of the individual emitters:
counts     = zeros(length(power),nEMITTERS);
satPower   = zeros(nEMITTERS,1);

%There are three simulation options above. If one is activated, here the
%corresponding parameters are changed such, that an array of values up to
%the maximum one set above is generated. If it is deactivated, the array
%will consist of the maximum value set above only. In that way you can
%decide if you want to probe only the value you choosed, or also values
%below that to check the influence on the result:
if SimPsatRanges
    pSatRangeMin = 0;
else
    pSatRangeMin = P_SAT_RANGE;
    SAT_RANGE_STEPS = 1;
end
if SimEmitterNumbers
    nMaxEmitters = nEMITTERS;
else
    nMaxEmitters = 1;
    nEMITTER_STEPS = 1;
end
if SimMaxPowers
    nMaxPowers = MAX_POWER_STEPS;
else
    nMaxPowers = 1;
end

%Finally initialize the array for the averages and standard deviations of
%the saturation powers obtained as the true values, from the standard
%saturation fit and from the logarithmic ensemble fit (defined below):
pSatRanges     = linspace(pSatRangeMin,P_SAT_RANGE,SAT_RANGE_STEPS);
nEmitterArray  = linspace(1,nMaxEmitters,nEMITTER_STEPS);
PsatTrueAvg    = zeros(nEMITTER_STEPS,length(pSatRanges),nMaxPowers);
PsatTrueStd    = zeros(nEMITTER_STEPS,length(pSatRanges),nMaxPowers);
PsatSatFitAvg  = zeros(nEMITTER_STEPS,length(pSatRanges),nMaxPowers);
PsatSatFitStd  = zeros(nEMITTER_STEPS,length(pSatRanges),nMaxPowers);
PsatLogFitAvg  = zeros(nEMITTER_STEPS,length(pSatRanges),nMaxPowers);
PsatLogFitStd  = zeros(nEMITTER_STEPS,length(pSatRanges),nMaxPowers);


%% ---- Perform simulation ------------------------------------------------

%Now iterate over all chosen simulation params and the number of runs
for iSatRange  = 1:length(pSatRanges)
    pSatRange  = pSatRanges(iSatRange);
    PsatTrue   = zeros(nRUNS,1);
    PsatSatFit = zeros(nRUNS,1);
    PsatLogFit = zeros(nRUNS,1);
    for nEmitterIndex = 1:length(nEmitterArray)
        nEmitter = round(nEmitterArray(nEmitterIndex));
        for iRun = 1:nRUNS
            %Create a saturation curve for each emitter, where its
            %saturation power is taken from a uniform distribution with the
            %parameters defined above. The initial function for the count
            %rate is obtained from a physical model where a uniform emitter
            %area is probed with a laser and has logarithmic behaviour:
            for iEmitter = 1:nEmitter
                satPower(iEmitter) = (rand*pSatRange) + P_SAT_MEAN - ...
                                      pSatRange/2;
                counts(:,iEmitter) = 1/log(2).*...
                                     log(1+power./satPower(iEmitter));
            end
            %Get now the mean to find the ensemble result:
            PsatTrue(iRun) = mean(satPower(satPower>0))/P_SAT_MEAN;
            meanCounts     = mean(counts,2);
            %Next define the two fit functions, where the first corresponds
            %to a single emitter, and the second to a uniformly distributed
            %emitters ensemble with saturation powers in a range of 2*Psat
            %around a mean value Psat:
            normSatFct = @(I_0,P_0,x) I_0*x./(x+P_0);
            logSatFct  = @(I_0,P_0,x) I_0/log(2)*(x.*log(1+2*P_0./x) + ...
                                      2*P_0.*log(1+x./P_0./2));
            %Set the fit bounds:
            lb = [0,0];
            ub = [2*P_SAT_MEAN,2*P_SAT_MEAN];
            %And perform the fits for each maximum power:
            for iMaxPower = 1:nMaxPowers
                maxPower = MAX_POWER/nMaxPowers*iMaxPower;
                [SatFit] = fit(power(power<maxPower)',...
                               meanCounts(power<maxPower),normSatFct,...
                               'StartPoint',[P_SAT_MEAN,P_SAT_MEAN],...
                               'Lower',lb,'Upper',ub); 
                [LogFit] = fit(power(power<maxPower)',...
                               meanCounts(power<maxPower),logSatFct,...
                               'StartPoint',[P_SAT_MEAN,P_SAT_MEAN],...
                               'Lower',lb,'Upper',ub); 
                %Extract the results:
                resSat = [SatFit.I_0,SatFit.P_0];
                resLog = [LogFit.I_0,LogFit.P_0];
                satCurve = normSatFct(SatFit.I_0,SatFit.P_0,power);
                logCurve =  logSatFct(LogFit.I_0,LogFit.P_0,power);
                PsatSatFit(iRun,iMaxPower) = SatFit.P_0/P_SAT_MEAN;
                PsatLogFit(iRun,iMaxPower) = LogFit.P_0/P_SAT_MEAN;
            end
            %Print here the current calculation to track the progress:
            Current_Calculation = ...
            strcat(sprintf('Psat range: %.0f/%.0f --- ',...
                            iSatRange, SAT_RANGE_STEPS),...
            sprintf('emitter number: %.0f/%.0f --- ',...
                            nEmitterIndex, nEMITTER_STEPS),...
            sprintf('run: %.0f/%.0f',iRun, nRUNS)) %#ok<NOPTS>
            %And potentially plot a sat curve, if the option is activated:
            if plotSatCurves
                figure
                hold on
                plot(power,meanCounts,'kx')
                plot(power,satCurve,'r--')
                plot(power,logCurve,'b--')
                title(sprintf('Saturation curve %.0f for %.0f emitters',...
                              iRun,nEmitter))
                xlabel('power / mW')
                ylabel('counts (normalized)')
            end
        end
        %Get the mean values of the parameters for all runs:
        PsatTrueAvg  (nEmitterIndex,iSatRange,:) = mean(PsatTrue);
        PsatSatFitAvg(nEmitterIndex,iSatRange,:) = mean(PsatSatFit);
        PsatLogFitAvg(nEmitterIndex,iSatRange,:) = mean(PsatLogFit);
        %And the corresponding standard deviations:
        PsatTrueStd  (nEmitterIndex,iSatRange,:)  = std(PsatTrue);
        PsatSatFitStd(nEmitterIndex,iSatRange,:)  = std(PsatSatFit);
        PsatLogFitStd(nEmitterIndex,iSatRange,:)  = std(PsatLogFit);
    end
end


%% ---- Plot results ------------------------------------------------------

%Finally plot the results for the three different simulation modes you have
%chosen. Always the normalized (to the true mean value defined in the 
%beginning) Psat values obtained are plotted against the parameter that is
%varied in the given simulation mode:

if SimPsatRanges
    figure
    errorbar(pSatRanges,squeeze(PsatTrueAvg(end,:,end)),...
        squeeze(PsatTrueStd(end,:,end))/sqrt(nRUNS),'r--')
    hold on
    errorbar(pSatRanges,squeeze(PsatSatFitAvg(end,:,end)),...
        squeeze(PsatSatFitStd(end,:,end))/sqrt(nRUNS),'b--')
    errorbar(pSatRanges,squeeze(PsatLogFitAvg(end,:,end)),...
        squeeze(PsatLogFitStd(end,:,end))/sqrt(nRUNS),'k--')
    xlabel('range of P_{sat}')
    ylabel('effective P_{sat} / mean P_{sat}')
    legend('true value','sat fit','log fit')
end

if SimEmitterNumbers
    figure
    errorbar(nEmitterArray,squeeze(PsatTrueAvg(:,end,end)),...
        squeeze(PsatTrueStd(:,end,end))/sqrt(nRUNS),'r--')
    hold on
    errorbar(nEmitterArray,squeeze(PsatSatFitAvg(:,end,end)),...
        squeeze(PsatSatFitStd(:,end,end))/sqrt(nRUNS),'b--')
    errorbar(nEmitterArray,squeeze(PsatLogFitAvg(:,end,end)),...
        squeeze(PsatLogFitStd(:,end,end))/sqrt(nRUNS),'k--')
    xlabel('number of emitters')
    ylabel('effective P_{sat} / mean P_{sat}')
    legend('true value','sat fit','log fit')
end

if SimMaxPowers
    figure
    powerStepsNorm = MAX_POWER/nMaxPowers/P_SAT_MEAN;
    maxPowerArray = powerStepsNorm:powerStepsNorm:MAX_POWER/P_SAT_MEAN;
    errorbar(maxPowerArray,squeeze(PsatTrueAvg(end,end,:)),...
        squeeze(PsatTrueStd(end,end,:))/sqrt(nRUNS),'r--')
    hold on
    errorbar(maxPowerArray,squeeze(PsatSatFitAvg(end,end,:)),...
        squeeze(PsatSatFitStd(end,end,:))/sqrt(nRUNS),'b--')
    errorbar(maxPowerArray,squeeze(PsatLogFitAvg(end,end,:)),...
        squeeze(PsatLogFitStd(end,end,:))/sqrt(nRUNS),'k--')
    xlabel('max Power simulated / mean P_{sat}')
    ylabel('effective P_{sat} / mean P_{sat}')
    legend('true value','sat fit','log fit')
end


%% ---- Protocol of updates -----------------------------------------------

% 25.01.21 (V1.1):  - Sections introduced
%                   - Structure improved
%
% 26.01.21 (V1.2):  - Option for changing maximum power introduced
%                   - "Current calculation" print out added
%                   - Appearance of the fits changed
%                   - Code commented
%
% 01.02.21 (V1.3):  - Warnings removed
