%
%,------------------------------------------------------------------------,
%| Simulation of FWHM of a scanned spot           (Lukas @27/01/21)  V1.2 |
%'------------------------------------------------------------------------'
%
%In this script the FWHM of a scanned emitter spot is calculated, as it
%appears from the scan. The laser is sampled over the spot and the
%resulting signal is collected. The uniformely implanted spot is then
%mapped with this scan. Inhere you can change the laser focus size and the
%emitter spot size, to see how this effects the result. Keep in mind, that
%in case you collect with the same objective as you excite, the focal spot
%width that you give in the following as w_0 relates to the width of the
%total PSF, which is the PSF of the laser times that of the emitter.

close all
clear variables


%% ---- Parameters to be predefined ---------------------------------------

r_0 = 500;      %Radius of the uniformely implanted spot
X_MAX = 2000;   %Maximum distance of Gauss from center that is considered
W_MAX = 500;    %Maximum Gaussian width considered
STEP  = 100;    %Steps in which Gaussian width is changed

%! If you choose a ratio W_MAX/STEP>10 the simulation might take long !


%% ---- Calculation of the appearing spot in the scan ---------------------

%The following function is the analytic integral in y-direction, of a 
%2D-Gaussian function over a circle, that is shifted along the x-axis. You
%see basically the error functions as the integrated Gauss function along
%one axis with the edges of the circle as argument (integral bounds). These
%are multiplied by the remaining part of the 2D Gauss, that is the x
%dependend one. You will need to repead the integration for the other axis,
%but since the bounds are dependend on x via the geometry of the circle,
%the integral becomes impossible to solve analytically. That is why the
%second integration will be done numerically in the following.

Fct = @(x_0,r_0,w_0,x) (erf(sqrt(r_0^2-(x-x_0).^2)/(sqrt(2)*w_0))-...
                        erf(-sqrt(r_0^2-(x-x_0).^2)/(sqrt(2)*w_0))).*...
                        exp(-x.^2/(2*w_0^2));
                    
%Initialize relevant arrays:
w_0     = zeros(1,W_MAX/STEP);          %Gaussian focus width
NumInt  = zeros(X_MAX+1,1);             %result of numeric integral
Spot    = zeros(X_MAX+1,W_MAX/STEP);    %resulting spot on the sample
numFWHM = zeros(1,W_MAX/STEP);          %numeric FWHM of that spot

%Iterate now over all width of the Gauss and integrate the function for
%each x-position over the circle, which is equivalent to convoluting the
%2D-Gauss with a circle along the x-axis. Calculate then the numerically
%derived FWHM of the resulting distribution.

for iGaussWidth=1:W_MAX/STEP
    w_0(iGaussWidth) = iGaussWidth*STEP;
    for x_0=0:X_MAX
        NumInt(x_0+1) = integral(@(x) Fct(x_0,r_0,w_0(iGaussWidth),x),...
                                          x_0-r_0,x_0+r_0);
    end
    Spot(:,iGaussWidth)  = NumInt;
    numFWHM(iGaussWidth) = 2*length(NumInt(NumInt>max(NumInt/2)));
end
  
  
%% ---- Plot the results --------------------------------------------------

figure('Position', [420 120 920 820])

subplot(2,1,1)
for iGaussWidth = 1:length(Spot(1,:))
    plot(Spot(:,iGaussWidth)/max(Spot(:,iGaussWidth)))
    hold on
end
ylim([0,1.1])
title('Scanned spot for different Gaussian width')
xlabel('radial distance from center of implantation / nm')
ylabel('normalized intensity')
legend(sprintf('%.0f nm',STEP),sprintf('%.0f nm',2*STEP),...
       sprintf('%.0f nm',3*STEP),'...')
   
subplot(2,1,2)
plot(w_0,numFWHM)
title('Dependence of spot FWHM on Gaussian width')
xlabel('Gaussian stdev / nm')
ylabel('spot FWHM')


%% ---- Protocol of updates -----------------------------------------------

% 27.01.21 (V1.2):  - Warnings removed.