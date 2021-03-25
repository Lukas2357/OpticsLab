%
%,------------------------------------------------------------------------,
%| Simulation of emitter spot scanning            (Lukas @27/01/21)  V2.1 |
%'------------------------------------------------------------------------'
%
%In this script a confocal scan of a spot of emitters is simulated. The
%emitters are arranged homogeniously in a given sized circle or Gaussian 
%distributed around a certain spot. A focused Gaussian beam by standard 
%objective is assumed and scanned over the spot. 

%The program will simulate the resulting intensity profile of the collected
%emission. This is done as follows: In the SimTotalPSF.m script, the total
%PSF (excitation * emission) of a single emitter is generated, including
%saturation effects. It is then fitted by a Gaussian profile and the
%parameters are saved for different saturation parameters. These are then
%loaded by this program. In addition a circular heaviside function or a 2D 
%radial Gaussian distribution is generated and then manually convoluted 
%with the Gauss. In that way the scanning profile of the spot is simulated.

close all
clear variables


%% ---- Optional functions ------------------------------------------------

showImage    = false;         %Show resulting image of the spots from scan
manualParams = true;          %Set focus params manually instead of loading
homogSpot    = false;         %Homogeneous emitter spot (else Gauss distr.)

%% ---- Define parameters -------------------------------------------------

RHO_MAX     = 1000;     %Max simulated radius from center of impl. (nm)
%Choose a radius, which covers the emitter and laser spot completely!
STEP        = 20;       %Steps in which spot is simulated (nm)
%A step value more than 100x less than RHO_MAX can take very long sim time!
SPOT_RADIUS = 500;      %Radius of implanted spot in homogeneous case (nm)
%This option is only relevant if you use homogeneous emitter distribution!
SIGMA_X     = 150;      %X-Stddev of Gaussian distributed emitter spot (nm)
SIGMA_Y     = 150;      %Y-Stddev of Gaussian distributed emitter spot (nm)
%These two are used if homogSpot is deactive. Attention: In order to allow
%the numerical fit fct, to properly describe the result associated with
%this Gaussian distribution, you need to adjust local params in that fct
%according to the values provided here, at the end of the script!

%If you activate manualParams you can provide them here:
if manualParams
    %The SatParam defines how strong the emitters are driven. If it is 1,
    %the peak intensity in the laser focus matches the saturation power of
    %the emitter, for >1 emitters will saturate, affecting the spot size:
    SatParam = linspace(0,10,30);
    %The PeakStd gives the standard deviation of the lateral focus size. It
    %is typically defined by the objectives NA and filling and will
    %strongly affect the appearance of the emitter spot in the scan:
    PeakStd  = linspace(150,150,30);
else
    cd Data %#ok<*UNRCH>
    load('SatParam')        %Loads saturation parameter from data folder
    SatParam = SatParam';   %Reshape SatParam array to match it with rest
    load('PeakStd')         %Loads stddev of focal spot from data folder
    cd ..
end


%% ---- Initialize variables ----------------------------------------------

rho   = -RHO_MAX:STEP:RHO_MAX;             %Simulated radial intervall
rho_0 =  RHO_MAX/STEP+1;                   %Center of the intervall (rho=0)
mesh  = zeros(length(rho),length(rho),2);  %Initial. sample positions mash

for i = 1:length(rho)
    mesh(i,:,1) = rho;      %Use the rho array to define a cartesian cs
    mesh(:,i,2) = rho;
end
x = mesh(:,:,1);            %And extract the corresponding x and y matrices
y = mesh(:,:,2);

%Spatial intensity distribution (2D) of the focus in each space point (2D), 
%resulting in a 4D object called FocalSpot:
FocalSpot = zeros(length(x),length(y),2*RHO_MAX/STEP+1,2*RHO_MAX/STEP+1);
%Spatial emitter distribution (2D) in a fixed position, define as a
%homogeneous circle of emitters via heaviside function (2D object) or via a
%two dimensional Gaussian distribution:
if homogSpot
    EmitterSpot = (1-heaviside(sqrt(x.^2+y.^2)-SPOT_RADIUS));
else
    EmitterSpot = exp(-(x/SIGMA_X).^2/2-(y/SIGMA_Y).^2/2);
end
%Multiplying FocalSpot with EmitterSpot for each position results in a 2D
%distribution of intensity for each point, that is in total a 4D object
%called Counts in here:
Counts = zeros(length(x),length(y),length(PeakStd));
%To construct the Image of the spot from this, sum up the distribution in
%each position, the result is a 2D object called Image here. It will be
%created for each satParam, so the overall object is 3D. The same is done
%for the array used for fitting later (FitImage):
Image    = zeros(2*RHO_MAX/STEP+1,2*RHO_MAX/STEP+1,length(PeakStd));
FitImage = zeros(2*RHO_MAX/STEP+1,2*RHO_MAX/STEP+1,length(PeakStd));
%For each saturation parameter the maximum count rate in the image is
%identified and saved in MaxCounts:
MaxCounts = 1:length(PeakStd);


%% ---- Calculate convolution of focal spot with emitter spot -------------

%Now iterate over all elements of PeakAmp or PeakStd, respectively, to
%generate results for the different saturation paramters:
for iAmp = 1:length(PeakStd)
    %Then iterate over all matrix elements in the x- and y-position matrix
    %to generate the effective count rate in each position. This is done as
    %follows: The emitter circle is fixed and the position of the focal
    %spot is varied over all positions in the chosen area (up to RHO_MAX).
    %In each position the product of the circle and the Focus is calculated
    %and the resulting profile is integrated to find the effective counts
    %in this position. Afterwards all these values for all positions are
    %put together in a matrix to find the resulting image.
    for j = 1:2*RHO_MAX/STEP+1
        for k = 1:2*RHO_MAX/STEP+1
            x_0 = (j-1)*STEP-RHO_MAX;
            y_0 = (k-1)*STEP-RHO_MAX;
            FocalSpot = 1./(1+exp(((x(1,:) - x_0)/PeakStd(iAmp)).^2/2 +...
                                  ((y(:,1) - y_0)/PeakStd(iAmp)).^2/2)/...
                                  SatParam(iAmp));
            Counts = EmitterSpot.*FocalSpot;
            if x_0==0 && y_0==0
                centralFocalSpot = FocalSpot;
                centralCounts = Counts;
            end
            Image(j,k,iAmp) = sum(Counts(:));
        end
    end
    MaxCounts(iAmp) = max(max(Image(:,:,iAmp)));


    %% ---- Fit image spot and plot result --------------------------------

    %---- Standard fitting routine using 2D Gauss function ----------------
    GuessParams = [1,SPOT_RADIUS,SPOT_RADIUS];
    [Fit_result,resnorm,residual,exitflag,output,lambda,jacobian]  = ...
           lsqcurvefit(@D2GaussFunction,GuessParams,mesh,Image(:,:,iAmp)); 
    FitImage(:,:,iAmp) = D2GaussFunction(Fit_result,mesh);
    FWHM = 2*sqrt(log(2))*(Fit_result(2));  %Spatial FWHM of fitted Gauss
    %----------------------------------------------------------------------
    
    %Plot images of spot and laser focus in the scan, if activated:
    if showImage
        figure('Position', [420 120 1020 820])
        subplot(2,2,1)
        surf(x,y,centralFocalSpot)
        colormap(infernoColorMap)
        title  ('Focal spot')
        xlabel ('y-position / nm')
        ylabel ('x-position / nm')
        zlabel ('normalized intensity')
        subplot(2,2,2)
        surf(x,y,EmitterSpot)
        title  ('Emitter spot')
        xlabel ('y-position / nm')
        ylabel ('x-position / nm')
        zlabel ('normalized intensity')
        subplot(2,2,3)
        surf(x,y,centralCounts)
        title  ('Convolution at the center')
        xlabel ('y-position / nm')
        ylabel ('x-position / nm')
        zlabel ('normalized intensity')
        subplot(2,2,4)
        surf(x,y,Image(:,:,iAmp),'FaceAlpha',0.9)
        hold on
        surf(x,y,FitImage(:,:,iAmp),'LineStyle','-','FaceAlpha',0)
        title  ('Image and Gaussian fit')
        xlabel ('y-position / nm')
        ylabel ('x-position / nm')
        zlabel ('normalized intensity')
    end
    fprintf('Result of SatParam %.0f/%.0f calculated',iAmp,length(PeakStd))
end


%% ---- Investigate saturation of image spot ------------------------------

%The image spot shows a saturation behaviour according to the saturation of
%the emitters. But due to the convolution the saturation of the image spot
%could actualy have a different saturation parameter. This is investigated
%by fitting a saturation curve (SatFct) to the maximum counts of the images
%as well as a logarithmic fct and a numerically integrated fct. For more
%details on these models see Appendix D, Masterthesis Lukas Hunold.

SatFct = @(sat_amp,p_sat,x) sat_amp*x./(x+p_sat);
LogFct = @(log_amp,p_sat,x) log_amp*log(1+x./p_sat)/log(2);
lb = [0,0];
ub = [10*max(MaxCounts),10];
disp('Evaluating saturation fits...')
[SatResult] = fit(SatParam(SatParam>0)',MaxCounts(SatParam>0)',SatFct,...
                  'StartPoint',[max(MaxCounts),max(SatParam)],...
                  'Lower',lb,'Upper',ub);
disp('First done...')
[LogResult] = fit(SatParam(SatParam>0)',MaxCounts(SatParam>0)',LogFct,...
                  'StartPoint',[max(MaxCounts),max(SatParam)],...
                  'Lower',lb,'Upper',ub);
disp('Second done...')
[IntResult] = lsqcurvefit(@fitFunc,[9,1],SatParam(SatParam>0)',...
                           MaxCounts(SatParam>0)',lb,ub);
disp('Almost ready...')
%Generate the corresponding curves:
SatCurve = SatResult.sat_amp*SatParam./(SatParam+SatResult.p_sat);
LogCurve = LogResult.log_amp*log(1+SatParam./LogResult.p_sat)/log(2);
IntCurve = fitFunc(IntResult,SatParam);


%% ---- Plot the maximum counts of the images together with the fits ------

figure('Position',[681,652,467,327])
normFact = min(MaxCounts(SatParam>1));
plot(SatParam,MaxCounts/normFact,'k.','MarkerSize',10)
hold on
plot(SatParam,SatCurve/normFact,'-','LineWidth',1.5)
plot(SatParam,LogCurve/normFact,'--','LineWidth',1.5)
plot(SatParam,IntCurve/normFact,':','LineWidth',1.5)
title  ('Saturation of emitter ensemble')
xlabel ('saturation parameter (P/P_{sat})')
ylabel ('normalized central spot counts')
str = strcat('Initial    P_{sat}  = 1.00','\newline',...
      sprintf('SatFct  P_{sat}  = %.2f',SatResult.p_sat),'\newline',...
      sprintf('LogFct  P_{sat}  = %.2f',LogResult.p_sat),'\newline',...
      sprintf('NumFct P_{sat} = %.2f',IntResult(2)));
annotation('textbox',[0.58 0.45 0.07 0.02],'String',{str},...
           'FitBoxToText','on');
legend('spot peak counts','saturation fct fit','logarithm fct fit',...
       'numerical fct fit','Location','SouthWest')
set(gca,'FontSize',12)


%% ---- Used functions ----------------------------------------------------

%This fct models the spot as it appears in the scan:
function F = D2GaussFunction(x,xdata)
    F = x(1).*exp(-((xdata(:,:,1)).^2/(2*x(2)^2) + ...
                    (xdata(:,:,2)).^2/(2*x(3)^2)));
end

%This fct is a fit fct for the saturation of a Gaussian distributed
%emitter ensemble:
function fitData = fitFunc(x,xdata)
    fitData = zeros(size(xdata));
    %Enter here the stddev of laser spot (sigma) and Gaussian emitter
    %distribution in x and y (sigma_x, sigma_y). Notice that these values
    %also need to be known when this fct is used to analyse a measurement!
    sigma = 150;
    sigma_x = 150;
    sigma_y = 150;
    for i = 1:length(xdata)
        fitData(i) = centerConvolution(x,xdata(i),sigma,sigma_x,sigma_y);
    end
end

%This fct is the underlaying numerical model for the above fit fct:
function int_res = centerConvolution(x,I_0,sigma,sigma_x,sigma_y)
    amp     = x(1);
    I_sat   = x(2);
    spot    = @(r,phi) amp.*exp(-r.^2.*(cos(phi).^2./2./sigma_x.^2 + ...
                                      sin(phi).^2./2./sigma_y.^2));
    integr  = @(r,phi) sigma.^-2.*r.*spot(r,phi).*...
                     (1+I_sat./I_0.*exp(r.^2./sigma.^2/2)).^-1;
    int_res = integral2(integr,0,5000,0,2*pi);
end


%% ---- Protocol of updates -----------------------------------------------

% 26.08.20 (V1.5):  - Sections introduced.
% 28.01.21 (V2.1):  - Documentation and structure improved
%                   - Manual params option introduced
%                   - Amplitude parameter removed (obsolete)
%                   - Log and Numerical fit functions introduced
%                   - Gaussian spot distribution introduced
%                   - Appearance of the plots optimized
