%
%,------------------------------------------------------------------------,
%| Plot and Fit axial intensity scans     (Lukas Hunold @17/01/21)  V 1.1 |
%'------------------------------------------------------------------------'
%
%In this script an intensity scan in axial (z) direction can be analysed. 
%You can do different z-scans with the APDs and the corresponding LabView
%program and save the data in a txt file. If you measure the reflection of
%a laser, as well as the signal of the emitters in that way, you will get
%peaks, that can be fitted with Gaussian or Lorentzian functions. 

close all
clear variables

%Here first a Lorentz fit is done as rough estimation. To take into account
%the antiysmmetric shape then, a Gaussians fit is performed for left and
%right part of each peak. With this you can extract the parameters and find
%out more about the size of your laser spot in z-direction as well as the
%distribution of your emitter and their depth below the surface. Save the 
%measurements as "measurement_name".txt and "measurement_name"-laser.txt, 
%where you set the variable "measurement_name" here:
measurement_name = 'test';
%The script will also work if you measured each point in z multiple times,
%it will then just take them, calculate the average and use it for the fit.
%Decide now if you want to save the results and then run the script:
SAVE_RESULT = true;


%% ---- Data recording and ordering ---------------------------------------

%Take data from txt files (columns of z and counts). The laser file
%corresponds to surface reflection and the other to the recorded signal:
cd Data
Data_Laser = importdata(strcat(measurement_name, '-laser.txt'));
Data_Signal  = importdata(strcat(measurement_name, '.txt'));
cd ..

%Extract z and counts separately for laser reflection and signal:
z_axis = Data_Laser(:,1);
Laser = Data_Laser(:,2)/max(Data_Laser(:,2));
Signal = Data_Signal(:,2)/max(Data_Signal(:,2));

%Initialize the averaged vectors:
for z_index = 1:length(z_axis)
    if z_axis(z_index+1) == z_axis(z_index)
        continue
    else
        break
    end
end
nPonits_per_z = z_index;
z_axis        = unique(z_axis);
Signal_mean   = 1:length(z_axis);
Laser_mean    = 1:length(z_axis);
Signal_std    = 1:length(z_axis);
Laser_std     = 1:length(z_axis);

%Peform the averaging over certain range of values:
for i = 1:length(Signal)/nPonits_per_z
    Signal_mean(i) = mean(Signal(nPonits_per_z*(i-1)+1:i*nPonits_per_z));
    Signal_std(i) = std(Signal(nPonits_per_z*(i-1)+1:i*nPonits_per_z));
end
for i = 1:length(Signal)/nPonits_per_z
    Laser_mean(i) = mean(Laser(nPonits_per_z*(i-1)+1:i*nPonits_per_z));
    Laser_std(i) = std(Laser(nPonits_per_z*(i-1)+1:i*nPonits_per_z));
end
Signal_mean = Signal_mean';
Laser_mean = Laser_mean'; 

%Calculate the difference of the maxima in z for comparison:
difference = z_axis(Laser_mean == max(Laser_mean))- ... 
             z_axis(Signal_mean == max(Signal_mean));

%And get the maximum values to have peak position and hight:
[LaserMax,iLaserMax] = max(Laser_mean);
[SignalMax,iSignalMax] = max(Signal_mean);

%In addition define left and right flanks of laser and emitter peak:
lz1 = z_axis(1:iLaserMax);
rz1 = z_axis(iLaserMax:end);
lz2 = z_axis(1:iSignalMax);
rz2 = z_axis(iSignalMax:end);
lLaser = Laser_mean(1:iLaserMax);
rLaser = Laser_mean(iLaserMax:end);
lSignal = Signal_mean(1:iSignalMax);
rSignal = Signal_mean(iSignalMax:end);


%% ---- Performing the Lorentz fit ----------------------------------------

%Define Lorentz fit function:
Lorentz = @(mvLor,ampLor,widthLor,x) ...
          ampLor .* widthLor^2./((x-mvLor).^2+widthLor^2);
    
%Perform the fit and extract the parameters:
[LorRes,LorGof,LorInfoout] = fit(z_axis,Signal_mean,Lorentz,...
                    'StartPoint',[mean(z_axis),1,max(z_axis)-min(z_axis)]);
LorentzParams  = coeffvalues(LorRes);
LorentzConf    = confint(LorRes);
LorentzWLc     = LorRes.mvLor;
LorentzWLcStd  = LorRes.mvLor - LorentzConf(1,1);
LorentzFWHM    = LorRes.widthLor*2;
LorentzFWHMStd = 2*(LorRes.widthLor - LorentzConf(1,3));
LorentzProfile = LorRes.ampLor*LorRes.widthLor^2 ./ ...
                 ((z_axis-LorRes.mvLor).^2+LorRes.widthLor^2);

%Perform the fit and extract the parameters:
[LorRes1,LorGof1,LorInfoout1] = fit(z_axis,Laser_mean,Lorentz,...
                    'StartPoint',[mean(z_axis),1,max(z_axis)-min(z_axis)]);
LorentzParams1  = coeffvalues(LorRes1);
LorentzConf1    = confint(LorRes1);
LorentzWLc1     = LorRes1.mvLor;
LorentzWLcStd1  = LorRes1.mvLor - LorentzConf1(1,1);
LorentzFWHM1    = LorRes1.widthLor*2;
LorentzFWHMStd1 = 2*(LorRes1.widthLor - LorentzConf1(1,3));
LorentzProfile1 = LorRes1.ampLor*LorRes1.widthLor^2 ./ ...
                 ((z_axis-LorRes1.mvLor).^2+LorRes1.widthLor^2);

             
%% ---- Plotting the Lorentz results --------------------------------------

%Plot the spectrum with the fit:
figure('Position',[1050,550,800,400])
errorbar(z_axis,Laser_mean,Laser_std/sqrt(nPonits_per_z),...
         '.k','LineWidth',1.5)
hold on
errorbar(z_axis,Signal_mean,Signal_std/sqrt(nPonits_per_z),...
         '.b','LineWidth',1.5)
plot(z_axis,LorentzProfile,'-r','LineWidth',1.5)
plot(z_axis,LorentzProfile1,'-r','LineWidth',1.5)
title ('Z scans of laser and emitter signal with Lorentz fit')
xlabel ('z-position / \mum')
ylabel ('normalized counts')
xlim   ([min(z_axis),max(z_axis)])
ylim   ([0,1.3])
legend('off')
dim = [0.15 0.6 0.3 0.3];
str = {strcat('Laser:','\newline','z_0 = ', ...
       sprintf(' (%.2f +/- %.2f)', ...
       LorentzWLc,LorentzWLcStd)), ...
       sprintf('FWHM = (%.2f +/- %.2f)', ...
       abs(LorentzFWHM), abs(LorentzFWHMStd))};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
dim1 = [0.68 0.6 0.3 0.3];
str = {strcat('Signal:', '\newline', 'z_0 = ', ...
       sprintf(' (%.2f +/- %.2f)', ...
       LorentzWLc1,LorentzWLcStd1)), ...
       sprintf('FWHM = (%.2f +/- %.2f)', ...
       abs(LorentzFWHM1),abs(LorentzFWHMStd1))};
annotation('textbox',dim1,'String',str,'FitBoxToText','on');
set(gca,'FontSize',12)


%% ---- Save Lorentz result -----------------------------------------------

cd 'Result'
if SAVE_RESULT
    saveas(gcf, strcat(measurement_name, '-zscanLorentzFit'), 'png');
end
cd ..


%% ---- Perform the two sides Gauss fits ----------------------------------

%Now use the left and right flanks of the peaks defined before to do the
%antisymmetric Gaussian fits:

widthUb = max(z_axis)-min(z_axis);

FitGauss = lsqcurvefit(@Gauss,[1,LorentzWLc,1,0],lz1,lLaser, ...
                              [0,min(z_axis),0,0], ...
                              [LaserMax*1.1,max(z_axis),widthUb,1]); 
lLaser = Gauss(FitGauss,lz1);
Result(:,1) = FitGauss;

FitGauss = lsqcurvefit(@Gauss,[1,LorentzWLc,1,0],lz2,lSignal, ...
                              [0,min(z_axis),0,0], ...
                              [SignalMax*1.1,max(z_axis),widthUb,1]);
lSignal = Gauss(FitGauss,lz2);
Result(:,2) = FitGauss;

FitGauss = lsqcurvefit(@Gauss,[1,LorentzWLc,1,0],rz1,rLaser, ...
                              [0,min(z_axis),0,0], ...
                              [LaserMax*1.1,max(z_axis),widthUb,1]);
rLaser = Gauss(FitGauss,rz1);
Result(:,3) = FitGauss;

FitGauss = lsqcurvefit(@Gauss,[1,LorentzWLc,1,0],rz2,rSignal, ...
                              [0,min(z_axis),0,0], ...
                              [SignalMax*1.1,max(z_axis),widthUb,1]);
rSignal = Gauss(FitGauss,rz2);
Result(:,4) = FitGauss;


%% ---- Plot two sided Gauss results --------------------------------------

%And finally plot the results of these fits:

figure('Position',[1050,50,800,400])
errorbar(z_axis,Laser_mean,Laser_std/sqrt(nPonits_per_z),...
         '.k','LineWidth',1.5)
hold on
errorbar(z_axis,Signal_mean,Signal_std/sqrt(nPonits_per_z),...
         '.b','LineWidth',1.5)
hold on
plot(lz1,lLaser,'-.r','LineWidth',2)
hold on
plot(rz1,rLaser,'--r','LineWidth',2)
hold on
plot(lz2,lSignal,'-.r','color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
plot(rz2,rSignal,'--r','color',[0.6350 0.0780 0.1840],'LineWidth',2)
title ('Z scans of laser and emitter signal with two sided Gaussian fit')
xlabel('depth / \mum')
ylabel('normalized intensity')
legend('Reflected laser','Emitter signal', 'left Laser fit', ...
       'right Laser fit','left Signal fit','right Signal fit')
set(gca,'FontSize',12)


%% ---- Save two sided Gauss results --------------------------------------

cd 'Result'
if SAVE_RESULT
    saveas(gcf, strcat(measurement_name, '-zscanGaussFit'), 'png');
end
cd ..


%% ---- Used functions ----------------------------------------------------

function F = Gauss(x,xdata)
 F = x(1).*exp(-((xdata-x(2))/(2*x(3))).^2)+x(4);
end


%% ---- Protocol of updates -----------------------------------------------

% 17.01.21 (V1.1):  - Sections introduced
%                   - Two sided Gauss fits added
