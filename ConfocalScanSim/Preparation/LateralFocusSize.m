%
%,------------------------------------------------------------------------,
%| Calculating lateral width of focal spots       (Lukas @27/01/21)  V1.4 |
%'------------------------------------------------------------------------'
%
%In this script a focused Gaussian beam is simulated for different 
%objective fillings (f0). Field and intensity distributions are calculated
%around the focus point and the lateral width of that focus is determined.

close all
clear variables


%% ---- Define parameters -------------------------------------------------

%First set here the parameters you use in the experiment:
NA    = 0.95;                %Objective numerical aperture
WL    = 656;                 %Excitation wavelength
P_INC = 8.13;                %Power before objective in mW
P_EFF = 3.31;                %Power after objective in mW
TRANS = 0.9;                 %Fraction of transmitted power in objective

%The corresponding objective filling is calculated out of that:
real_f0 = (log(P_INC/(P_INC-P_EFF/TRANS))/2)^(-1/2);

%For comparison also other fillings are simulated in the following, here
%you can choose which:
MIN_f0_SIM   = 0.5;     %minimum simulated filling factor
MAX_f0_SIM   = 3;       %minimum simulated filling factor
f0_STEPS_SIM = 20;      %number of steps simulated between the two

%Finally set the parameters of the simulation to get a suitable compromise 
%between precision and time consumption:
RHO_MAX  = 1500;              %Maximum distance for calculating PSF in nm
STEP_RHO = 2;                 %Steps for calculating PSF in nm


%% ---- Calculate spatial distribution of the intensity -------------------

%First use the real filling, the other ones are automatically added below.
%Also initialize the field and parameter arrays:
f0       = real_f0;                        
I_00     = zeros(2*RHO_MAX/STEP_RHO+1,1);
PeakAmp  = zeros(f0_STEPS_SIM,1);
PeakStd  = zeros(f0_STEPS_SIM,1);
PeakFWHM = zeros(f0_STEPS_SIM,1);

%Open figure to plot the field distribution around the focus:
figure('Position', [1250 550 650 430])
%And start the calculation:
for f0Index=1:f0_STEPS_SIM
    for iRho=1:2*RHO_MAX/STEP_RHO+1
        radius = (iRho-1)*STEP_RHO-RHO_MAX;     
        %This is the fct under the integral depending on the NA, the
        %objective filling, the WL and the radius, leading to the field:
        f   = @(x) (1-exp(-4./(2*f0.^2))).*NA^2*exp(-x.^2/f0.^2).*...
                    x.*((1-NA^2*x.^2).^(1/4)+(1-NA^2*x.^2).^(-1/4)).*...
                    besselj(0,2*pi/WL*NA*x*radius);
        %I_00 is the integral over this fct with respect to x, from 0 to 1:
        I_00(iRho) = integral(f,0,1);
    end
    %You can now calculate the focus intensity by normalizing and squaring:
    normInt = (I_00/max(I_00)).^2;


    %% ---- Fitting the PSF and extracting the parameters -----------------

    %Define a Gaussian fit function
    Gauss = @(amp_G,std,x) amp_G.*exp(-((x)/std).^2/2);
    %as well as a vector for the radial position
    rho = -RHO_MAX:STEP_RHO:RHO_MAX;
    %and perform the spatial Gaussian fit of the total PSF
    [resultGaussFit,G_gof,G_infoout] = fit(rho',normInt,Gauss,...
                                           'StartPoint', [1,200]);
    %to be able to identify relative peak values and the width of curves:
    PeakAmp(f0Index) = resultGaussFit.amp_G;
    PeakStd(f0Index) = resultGaussFit.std;
    PeakFWHM(f0Index) = 2*sqrt(log(2)*2)*PeakStd(f0Index);


    %% ---- Plot the total PSF for different powers -----------------------

    plot(rho',I_00)
    title  ('Lateral focus field for different filling factors')
    xlabel ('position / nm')
    ylabel ('normalized intensity')
    xlim   ([-1500 1500])
    ylim   ([-0.2 1])
    hold on

    %Before restarting the loop, set the filling factor to the next step:
    f0 = MIN_f0_SIM + (MAX_f0_SIM - MIN_f0_SIM)/f0_STEPS_SIM*f0Index;  
end


%% ---- Plot final results ------------------------------------------------

%Generate an array for the filling factors plottet in the following:
f0_array = linspace(MIN_f0_SIM,MAX_f0_SIM,f0_STEPS_SIM); 
%And finally plot the focus FWHM depending on objective filling:
figure
plot(f0_array(2:end),PeakFWHM(2:end),'LineWidth',2)
hold on
plot([0,f0_array(end)],[PeakFWHM(1),PeakFWHM(1)],'LineWidth',2)
xlabel ('filling factor f_0')
ylabel ('lateral focus FWHM / nm')
xlim   ([min(f0_array),max(f0_array)])
ylim   ([0.9*min(PeakFWHM),max(PeakFWHM)])
set(gca,'FontSize',12)
legend('theoretical FWHM for NA = 0.95',...
       'FWHM for real case (f_0 = 1.82)')
str= ['P_{incident} = ',num2str(P_INC,3), ' mW \newline' , ...
      'P_{effective} = ',num2str(P_EFF,3), ' mW \newline \newline' , ...
      'Objective filling factor  = ',num2str(f0_array(1),3), '\newline', ...
      'Peak FWHM = ',sprintf('%.0f',roundn(PeakFWHM(1),1)),' nm'];
annot= annotation('textbox',[0.5 0.76 0.07 0.02], 'String',{ str},...
        'FitBoxToText','on','BackgroundColor','w');

    
    %% ---- Protocol of updates -----------------------------------------------

% 27.01.21 (V1.4):  - Warnings removed
%                   - Sections introduced
%                   - Appearance and documentation improved
