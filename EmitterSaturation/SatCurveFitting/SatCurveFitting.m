%
%,------------------------------------------------------------------------,
%| Analysis of saturation curves          (Lukas Hunold @01/02/21)  V 1.3 |
%'------------------------------------------------------------------------'
%
%In this script a saturation measurement of a (single) emitter is analysed. 
%Perform an axial scan for different excitation powers and record the
%emitter signal in a txt file. Please name the file according to the power 
%level that you used, e.g. for level 6.5 please choose "0650.txt" and for 
%level 10 it is "1000.txt". The program will find the maximum count rate
%in axial direction and take it as the emitter signal for that laser power.
%Then the curve for all powers is generated and fitted with the standard
%saturation function. A background correction is possible when providing a
%similar measurement for the background signal.

close all
clear variables


%% ---- Pre-settings ------------------------------------------------------

%Provide here the name of the folder within this folder, where the signal
%data txt files are stored:
Signal_Folder  = 'test_signal';
%If you want to also analyse and subtract a background, activate:
BgCorrection   = 1;
%And provide here the folder of the background files:
BG_Folder      = 'test_bg';
%Typically the uncertainties of the bg fit can be very large, since the 
%axial response is more or less flat. To get proper errors on that curve
%you should therefore estimate them beforehand and manually set them here:
REL_BG_ERROR   = 0.02;
%If you want to see plots of the axial scans, activate:
showAxialPlots = false;
%Finally the powers used for your measurement are needed. The lines below
%are an example of how to get the powers out of the levels, when the
%correlation between the two is known. You can do it as you wish, the
%important thing is that the array "Power" contains exactely those powers
%you used for the subsequent measurements:
STEPS = 0.25;
Level = 6.5:STEPS:10;
Power = (Level*1.4-9)/5*3.31;
%#ok<*UNRCH>


%% ---- Get counts for each power level -----------------------------------

%If a background correction is done, you generate two curves for the
%background and the signal, and combine both for the bg corrected signal:
if BgCorrection
    nCurves = 2;
else
    nCurves = 1;
end
%Initialize relevant arrays:
cd Bg
cd(BG_Folder)
txt_files = dir('**/*.txt');
PeakInt = zeros(length(txt_files),2);
PeakStd = zeros(length(txt_files),2);
cd ../..
%And start the data readout.
%If the correction is deactivated you skip the first curve of the bg:
for iCurve = (3-nCurves):2
    
    %% ---- Data readout and ordering -------------------------------------
    if iCurve == 1
        cd Bg
        cd(BG_Folder)
    else 
        cd Signal
        cd(Signal_Folder)
    end
    %Define struct array for all txt files in the folder:
    txt_files = dir('**/*.txt');
    %And iterate over all files:
    for iFile = 1:length(txt_files)
        %Load the z-axis and corresponding signal from each file:
        Data_Signal = importdata(txt_files(iFile).name);
        Data_Signal = Data_Signal.data;
        z_axis = Data_Signal(:,1);
        Signal = Data_Signal(:,2);
        %Check now if z values appear multiple times. This happens, if the
        %axial measurement takes multiple values at a given z position.
        %Typically the number of repetitions is the same for each z value,
        %but at the start of the measurement it might be different.
        %Therefore check the start of the measurement here:
        for z_index = 1:length(z_axis)
            if z_axis(z_index+1) == z_axis(z_index)
                continue
            else
                break
            end
        end
        %The loop is left as soon as the z values differ and the starting z
        %value is skipped in the following, since it causes problems when
        %not beeing repeated as many times as the other ones:
        z_axis = z_axis(z_index+1:end);
        Data_Signal = Data_Signal(z_index+1:end);
        %Now check the remaining z values:
        for z_index = 1:length(z_axis)
            if z_axis(z_index+1) == z_axis(z_index)
                continue
            else
                break
            end
        end
        %The loop is left when the z values differ. The number of
        %repetitions of a given z value is memorized here:
        nPointsZ = z_index;
        %Now the z array is made unique (duplicated entries removed):
        z_axis = unique(z_axis);
        %And an array for signal mean values at given z is initialized:
        Signal_mean = 1:length(z_axis);
        %The signal values for a given z can now be averaged and filled in:
        for i = 1:length(Signal)/nPointsZ
            Signal_mean(i) = mean(Signal(nPointsZ*(i-1)+1:i*nPointsZ));
        end
        %Finally the array is adjusted to allow fitting below:
        Signal_mean = Signal_mean';
        Signal_mean = Signal_mean(1:length(z_axis));

        %% ---- Performing the fit ----------------------------------------
        %The axial emitter signal typically follows a Gaussian fct:
        Gauss = @(mvGau,ampGau,stdGau,bgGau,x) ...
                        ampGau.*exp(-((x-mvGau)/(2*stdGau)).^2) + bgGau;
        %With this the obtained mean values are fitted:
        [GauRes,GauGof,GauInfoout] = fit(z_axis,Signal_mean,Gauss,...
           'StartPoint',[mean(z_axis),max(Signal_mean),...
           max(z_axis)-min(z_axis),0],'Lower',[min(z_axis),...
           0.5*max(Signal_mean),0,0],'Upper',[max(z_axis),...
           1.5*max(Signal_mean),max(z_axis)-min(z_axis),max(Signal_mean)]); 
        %And the relevant parameters are extracted:
        GaussParams  = coeffvalues(GauRes);
        GaussConf    = confint(GauRes,0.68);
        GaussAmp     = GauRes.ampGau;
        GaussAmpStd  = GauRes.ampGau - GaussConf(1,2);
        GaussFWHM    = 2*sqrt(log(2)*2)*GauRes.stdGau;
        GaussFWHMStd = (GauRes.stdGau-GaussConf(1,3))*...
                        GaussFWHM/GauRes.stdGau;
        GaussProfile = GauRes.ampGau.* ...
                       exp(-((z_axis-GauRes.mvGau)/...
                       (2*GauRes.stdGau)).^2) + GauRes.bgGau;
        %Important are the peak counts and standard deviation:
        if BgCorrection
            PeakInt(iFile,iCurve) = GaussAmp+GauRes.bgGau;
            PeakStd(iFile,iCurve) = GaussAmpStd;
        else
            PeakInt(iFile,1) = 0;
            PeakStd(iFile,1) = 0;
            PeakInt(iFile,iCurve) = GaussAmp+GauRes.bgGau;
            PeakStd(iFile,iCurve) = GaussAmpStd;
        end

        %% ---- Plotting the results --------------------------------------
        %Plot the spectrum with the fit:
        if showAxialPlots
            figure
            plot(z_axis,GaussProfile)
            hold on
            plot(z_axis,Signal_mean,'x')
            hold on 
            title ('Spot Autofocus')
            xlabel ('z-position / \mum')
            ylabel ('cps')
            xlim   ([min(z_axis),max(z_axis)])
            legend('off')
            dim = [0.55 0.6 0.3 0.3];
            str = {strcat(sprintf('Amp = (%.2f +/- %.2f) cps', ...
                   GaussAmp, GaussAmpStd))};
            annotation('textbox',dim,'String',str,'FitBoxToText','on');
            dim1 = [0.15 0.6 0.3 0.3];
        end
    end
    cd ../..
end


%% ---- Fit Saturation Curve ----------------------------------------------

%Initialize curve, fit result, parameter and fit data cell arrays:
Curve        = cell(3,1);
CurveStd     = cell(3,1);
satFitResult = cell(3,1);
P_sat        = cell(3,1);
I_sat        = cell(3,1);
errP_sat     = cell(3,1);
errI_sat     = cell(3,1);
satFit       = cell(3,1);

%Next the extracted values are used to generate and fit the sat curve:
for iCurve = (3-nCurves):3
    PeakStd(:,1) = PeakInt(:,1)*REL_BG_ERROR;
    %Now define the curves for bg, signal and bg corrected signal:
    if iCurve < 3
        Curve{iCurve} = PeakInt(:,iCurve);
        CurveStd{iCurve} = PeakStd(:,iCurve);
    else
        Curve{iCurve} = (PeakInt(:,2)-PeakInt(:,1));
        CurveStd{iCurve} = (PeakStd(:,2)+PeakStd(:,1));
    end
    %Define saturation fit function:
    SatFit = @(I_0,P_0,b,x) I_0*x./(x+P_0) + b;
    %Perform the fit:
    [satFitResult{iCurve}] = fit(Power',Curve{iCurve},SatFit, ...
        'StartPoint',[max(Curve{iCurve}),max(Power),min(Power)]); 
    %Extract the parameters:
    satFitConf          = confint(satFitResult{iCurve},0.68);
    P_sat{iCurve}       = satFitResult{iCurve}.P_0;
    I_sat{iCurve}       = satFitResult{iCurve}.I_0/2;
    errP_sat{iCurve}    = abs(P_sat{iCurve} - satFitConf(1,2));
    errI_sat{iCurve}    = abs(I_sat{iCurve} - satFitConf(1,1)/2);
    %Construct the fit curve:
    satFit{iCurve} = 2*I_sat{iCurve}*Power./(Power+P_sat{iCurve}) + ...
                     satFitResult{iCurve}.b;

end
%And plot the results:
figure('Position', [150 150 1050 650])
for iCurve = (3-nCurves):3  
    errorbar(Power,Curve{iCurve},CurveStd{iCurve},'LineStyle','none',...
             'Marker','.','MarkerSize',14, 'LineWidth',2)
    hold on
    plot(Power,satFit{iCurve},'LineWidth',2)
end
xlabel('power reaching the sample / mW')
ylabel('cps on one SPAD')
if BgCorrection
    legend('Background','Background Fit','Signal','Signal Fit',...
           'Bg corrected Signal','Background corrected Fit',...
           'Location','NorthWest')
else
    legend('Signal','Signal Fit','Bg corrected Signal',...
           'Background corrected Fit','Location','NorthWest')
end
set(gca,'FontSize',14)
dim = [0.142 0.67 0.4 0];
str = {sprintf('P_{sat} = (%.2f +/- %.2f) mW'  ,P_sat{3},errP_sat{3}), ...
       sprintf('I_{sat} = (%.0f +/- %.0f) cps',...
       roundn(I_sat{3},1),roundn(errI_sat{3},1))};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',14);


%% ---- Protocol of updates -----------------------------------------------

% 26.01.21 (V1.2):  - No background subtraction option added
%                   - Structure and Flexibility improved
%
% 01.02.21 (V1.3):  - Warnings removed
