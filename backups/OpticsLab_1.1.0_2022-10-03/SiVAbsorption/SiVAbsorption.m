%
%,------------------------------------------------------------------------,
%| SiV absorption cross section & emission (Lukas Hunold @30/01/21)  V1.3 |
%'------------------------------------------------------------------------'
%
%In this script the absorption of SiV centers is investigated with respect
%to the influence on effective emission count rates. The script can be
%reused for other emitters after slight modifications (see code below).

close all
clear variables


%% ---- Parameters to be predefined ---------------------------------------

%SI units used, except for lengths, these are given in nanometer always.

EXCITATION_WL = 656;
EMISSION_WL   = 738;

EMITTER_LIFETIME   = 1*10^-9;
EMISSION_LINEWIDTH = 5;

OBJ_INC_POW = 2*10^-3;          %Max laser power incident to the objective
OBJ_OUT_POW = 1*10^-3;          %Max laser power after the objective
OBJ_TRANSM  = 0.9;              %Fraction of power transmitted in the obj.
NA          = 0.95;             %NA of the objective

DIAMOND_TRANSM = 0.83;          %Fraction of power trasm. through surface

RHO_MAX = 2500;                 %Maximum focal radius simulated (nm)
STEP    = 1;                    %Steps to simulate radius (nm)


%% ---- Calculation of the peak intensity ---------------------------------

%First initialize radius array for focus calculation on the sample:
rho     = -RHO_MAX:STEP:RHO_MAX;
%Then an array for all effective excitation powers up to the maximum set:
excPower = DIAMOND_TRANSM*linspace(0,OBJ_OUT_POW,100);
%Calculate objective filling from incident and outgoing powers:
objFill = (log(OBJ_INC_POW/(OBJ_INC_POW-OBJ_OUT_POW/OBJ_TRANSM))/2)^(-1/2);

%Calculate the spatial extend in radial direction (radial symetric
%distribution) up to the maximum radius RHO_MAX in steps of STEP:
I_00 = zeros(2*RHO_MAX/STEP+1,1);
for i=1:2*RHO_MAX/STEP+1
    radius = (i-1)*STEP-RHO_MAX;     
    %This is the fct under the integral depending on the NA, the
    %objective filling, the WL and the radius:
    f_0   = @(z) (1-exp(-4./(2*objFill.^2))).*...
                NA^2*exp(-z.^2/objFill.^2).*...
                z.*((1-NA^2*z.^2).^(1/4)+(1-NA^2*z.^2).^(-1/4)).*...
                besselj(0,2*pi/EXCITATION_WL*NA*z*radius);
    %I_00 is the integral over this fct with respect to x, from 0 to 1:
    I_00(i) = integral(f_0,0,1);
end

%You can now calculate the PSF (spactial intensity distribution of the 
%focal spot) for r from 0 to RHO_MAX by squaring:
PSF = I_00.^2/max(I_00)^2;

%The following code gives you the Airy pattern minima and calculate
%from that the Airy disk (diameter):
minPSF = rho(PSF<10^-5);
for i = 1:length(minPSF)-1
    if minPSF(i) == minPSF(i+1)-STEP
        minPSF(i) = 0;
    end
end
minPSF    = nonzeros(minPSF);
AiryDisk  = min(minPSF(minPSF>0)) - max(minPSF(minPSF<0));
FocusFWHM = 2*max(rho(PSF>max(PSF)/2));

%Now integrate the PSF in polar coordinates by:
intPSF = sum(PSF'.*abs(rho))*pi;
%And calculate the intensity (in kW/cm^2 !!!) by normalizing the power to
%the integrated PSF and multiplying with the PSF itself:
intensity = excPower/intPSF.*PSF*10^11;  %10^11 for unit match
%The intensity at the center is now given (in kW/cm^2) by:
peakInt = max(intensity);
%If you would take a uniform distribution of the intensity over the airy
%disk, the intensity would be
uniformInt = excPower/(pi*(AiryDisk/2)^2)*10^11;
%Therefore you can calculate the peak intensity by using a uniform
%intensity and scaling afterwards by the factor:
intRatio = peakInt./uniformInt;

figure
plot(rho,intensity,'LineWidth',1.5)
title ('Intensity distributions of focal spot')
xlabel ('radius / nm')
ylabel ('intensity / (kW/cm^2)')
legend('off')
legend(sprintf('%.4f mW',excPower(1)*10^3), ...
       sprintf('%.4f mW',excPower(2)*10^3), ...
       sprintf('%.4f mW',excPower(3)*10^3),'...')
set(gca,'FontSize',12)
    

%% ---- Calculation of the absorption cross section -----------------------

%The absorption cross section is proportional to the ratio of broadened and
%natural linewidth, which are calculated as follows:
naturalLinewidth    = 1/(2*pi*EMITTER_LIFETIME);
broadenedLinewidth  = (3*10^8/(EMISSION_WL-EMISSION_LINEWIDTH/2)-...
                       3*10^8/(EMISSION_WL+EMISSION_LINEWIDTH/2))*10^9;
                   
%Additional factors to be taken into account are the following ones, which
%are from literature (you might change them, if you know better):
a_NR = 0.12;    %Fraction of non-radiative decay
a_DW = 0.7;     %Debay-Waller-Factor
a_OR = 1/6;     %Off resonance correction

%Absorption cross section (in cm^2 !!!) taking into account broadening and 
%factors above:
AbsCrossSection = 3*(EMISSION_WL*10^-9)^2/(2*pi) * ...
                  naturalLinewidth/broadenedLinewidth * ...
                  a_NR*a_DW*a_OR*10^4;
    
              
%% ---- Calculation of the absorbed photons -------------------------------

%Calculate the number of photons absorbed by the emitter placed in the
%center of the focus per second:
absorbablePhotons = AbsCrossSection*peakInt*1000/...
                    (2*pi*3*10^8/(EMISSION_WL*10^-9)*1.055*10^-34);
%The saturation parameter corresponding to that value is:
Saturation      = absorbablePhotons*EMITTER_LIFETIME;
%Which will result in a number of absorbed photons of:
absorbedPhotons = absorbablePhotons./(Saturation+1);

figure 
plot(excPower*1000,absorbablePhotons,'r-','MarkerSize',12,'LineWidth',2)
hold on
plot(excPower*1000,absorbedPhotons,'k.','MarkerSize',10,'LineWidth',1)
title ('Single emitter absorption')
xlabel('laser power after objective / mW')
ylabel('photons / cps')
xlim([min(excPower*1000),max(excPower*1000)])
legend('absorbable photons','absorbed photons','Location','NorthWest')
set(gca,'FontSize',12)


%% ---- Protocol of updates -----------------------------------------------

% 26.08.20 (V1.2):  - Sections introduced.
%
% 30.01.21 (V1.3):  - Warnings removed.
