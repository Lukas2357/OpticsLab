%
%,------------------------------------------------------------------------,
%| Simulate total PSF of an emitter               (Lukas @27/01/21)  V1.5 |
%'------------------------------------------------------------------------'
%
%In this script you can calculate and plot the total PSF of an emitter,
%that is excited and collected by the same objective, for different
%saturation parameters. The normalized peak intensity and the width of the
%total PSF are saved for each saturation parameter.

close all
clear variables


%% ---- Define parameters -------------------------------------------------

NA = 0.95;                  %Objective numerical aperture
WL = [656,738];             %Excitation and emission wavelength

P_SAT      = 10;            %Saturation power of the emitter in mW
MAX_POWER  = 25;            %Maximum power of the excitation in mW
POWER_STEP = 0.5;           %Steps of used powers for calculation in mW

RHO_MAX = 1000;             %Maximum distance for calculating PSF in nm
STEP    = 10;               %Steps for calculating PSF in nm


%% ---- Calculate objective filling ---------------------------------------

%To calculate the focal spot, the filling of the objective is important.
%It can be calculated from the incident power to the objective and the
%power afterwards by the formula below:

P_INC = 7.8;                %Power before objective in mW
P_EFF = 5.1;                %Power after objective in mW
TRANS = 0.9;                %Fraction of transmitted power in objective

objFill = (2*log(P_INC/(P_INC-P_EFF/TRANS)))^(-1/2);


%% ---- Calculate spatial distribution of the intensity -------------------

%We assume a Gaussian beam with wavelength WL(1), that is focused by an 
%objective with numerical aperture NA on an emitter. The emitter emits 
%light with wavelength WL(2) in a random direction, that is collected by 
%the same objective. To find the collected photons depending on the
%position of the objective relative to the emitter, one has to consider
%illumination and collection process. Therefore the number of emitted
%photons, that are collected by the objective in a certain position, is
%given by the product of excitation and emission point spread function.
%Effectively this is the square of the focal spot of the laser, but once
%with the wavelength of the emission. 
%In the following, the focal spot is described by the field distribution
%adapted from Novotny's Principles of Nano-Optics (p. 63). Neglecting
%higher order integrals, the intensity is proportional to the square of the
%I_00 integral given in the book. The function under this integral will be
%labeled as f in the following.

%Initialize arrays for the field (I_00) and the PSF:
I_00 = zeros(2*RHO_MAX/STEP+1,2);
PSF = zeros(2*RHO_MAX/STEP+1,2);

%Calculate now the focal spot for excitation and emission WL:
for wl = 1:2
    %Calculate the spatial extend in radial direction (radial symetric
    %distribution) up to the maximum radius RHO_MAX in steps of STEP:
    for i=1:2*RHO_MAX/STEP+1
        radius = (i-1)*STEP-RHO_MAX;     
        %This is the fct under the integral depending on the NA, the
        %objective filling, the WL and the radius:
        f   = @(x) (1-exp(-4./(2*objFill.^2))).*...
                    NA^2*exp(-x.^2/objFill.^2).*...
                    x.*((1-NA^2*x.^2).^(1/4)+(1-NA^2*x.^2).^(-1/4)).*...
                    besselj(0,2*pi/WL(wl)*NA*x*radius);
        %I_00 is the integral over this fct with respect to x, from 0 to 1:
        I_00(i,wl) = integral(f,0,1);
    end
    %You can now calculate the PSF (spactial intensity distribution of the 
    %focal spot) for r from 0 to RHO_MAX by normalizing and squaring:
    PSF(:,wl) = (I_00(:,wl)/max(I_00(:,wl))).^2;
end


%% ---- Saturation effects and total PSF ----------------------------------

%If the emitter is saturable by the laser intensity, the effective
%appearance of it in imaging or scanning will change. The reason for that
%is, that at high intensities the count rate will change more slowly with
%respect to the intensity. Therefore the change in the count rate with
%distance from emitter to objective will be smaller if the objective is
%close, since the intensity is high in that case. This effect only affects
%the excitation. It can be taken into account, by rescaling the excitation
%intensity such, that it represents the excitation of the emitter instead.
%That is easily done by inserting it in the saturation function of the
%emitter.

%Initialize the parameters (for faster computation)and figure:
nPowers  = (MAX_POWER/POWER_STEP);
PeakAmp  = zeros(nPowers,1);       %relative peak ampl. depending on power
PeakStd  = zeros(nPowers,1);       %relative peak width depending on power
SatParam = zeros(nPowers,1);       %satur. parameter depenfding on power
%And the power scale used:
Power   = POWER_STEP:POWER_STEP:MAX_POWER;

figure('Position', [1000 370 900 600])
%Calculate now for each power separately the total PSF:
for i = 1:nPowers
    %To do so notice first, that the saturation function is usually 
    %measured when the emitter is exactely in the focus. So the saturation 
    %power refers to the intensity in the center of the focus for that 
    %power. Now we want to construct the saturation function out of this, 
    %but for the emitter not exactely in the center of the focus. To do so,
    %we take the normal saturation curve with the used power and the 
    %sauration power and rescale the power by the relative intensities in 
    %the different distances from the center. This is valid, since the peak
    %intensity scales linear with the power. To find the total PSF we 
    %multiply with the emission PSF afterwards:
    totalPSF = PSF(:,1)*Power(i)./(PSF(:,1)*Power(i)+P_SAT).*PSF(:,2); 


    %% ---- Fitting the PSF and extracting the parameters ---------------------

    %We then define a Gaussian fit function
    Gauss = @(amp_G,std,x) amp_G.*exp(-((x)/std).^2/2);
    %as well as a vector for the radial position
    rho = -RHO_MAX:STEP:RHO_MAX;
    %and perform the spatial Gaussian fit of the total PSF
    [resultGaussFit,G_gof,G_infoout] = fit(rho',totalPSF,Gauss, ...
                                        'StartPoint', [1,200]);
    %to be able to identify relative peak values and width of the curves:
    PeakAmp(i) = resultGaussFit.amp_G;      %relative peak amplitude
    PeakStd(i) = resultGaussFit.std;        %peak width (stddev) in nm


    %% ---- Plot the total PSF for different powers ---------------------------

    plot(rho',totalPSF)
    %Calculate the saturation factor to display it in the figure:
    SatParam(i) = (i*POWER_STEP)/P_SAT;
    %Create an annotation for each power with its saturation factor:
    if i<9 %Show only first 8 elements for readability
        annotationHeight = 1.4*i/(MAX_POWER+1) + 0.2;
        str1   = ['S_0 = ',num2str(SatParam(i),3)];
        str2   = ['I_0 = ',num2str(PeakAmp(i),2)];
        str3   = ['\Delta I = ',num2str(PeakStd(i),4)];
        annot1 = annotation('textbox',[0.18 annotationHeight 0 0], ...
                            'String',{ str1}, 'FitBoxToText','on');
        annot2 = annotation('textbox',[0.6 annotationHeight 0 0], ...
                            'String',{ str2}, 'FitBoxToText','on');
        annot3 = annotation('textbox',[0.7 annotationHeight 0 0], ...
                            'String',{ str3}, 'FitBoxToText','on');         
    end
    title  ('Total PSF for different saturation parameters')
    xlabel ('position / nm')
    ylabel ('normalized intensity')
    xlim   ([-1000 1500])
    ylim   ([0 1])
    hold on

end


%% ---- Save relative amplitudes and width of the total PSFs --------------

cd ../Data
save('SatParam','SatParam')
save('PeakAmp','PeakAmp')
save('PeakStd','PeakStd')
cd ../Preparation


%% ---- Protocol of updates -----------------------------------------------

% 26.08.20 (V1.4):  - Sections introduced.
% 27.01.21 (V1.5):  - Warnings removed.

