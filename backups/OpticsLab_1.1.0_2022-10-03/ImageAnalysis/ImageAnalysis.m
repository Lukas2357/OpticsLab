%
%,------------------------------------------------------------------------,
%| Plot and Fit camera images              (Lukas Hunold @30/01/21)  V1.2 |
%'------------------------------------------------------------------------'
%
%In this script the image of a camera saved as .asc file is plottet and 
%eventually parts of it are fitted.

close all
clear variables


%% ---- Optional functions ------------------------------------------------

%First give the folder inside folder 'Data', where the files are located:
dataFolder = 'Test';
%And the name of the file:
fileName = 'test.asc';

%If you want to see the intensity map, activate:
generatePlots     = 1;
%If you want to see the full image in this case, activate:
generateFullImage = 0;
%Otherwise you have to give the x and y limits for the plot here:
PLOT_X_MIN = 150;
PLOT_X_MAX = 275;
PLOT_Y_MIN = 212;
PLOT_Y_MAX = 275;
%Here you can set the plot dimension for visualization:
PLOT_DIM = 3;   %options: 2 for 2D plot, 3 for 3D plot
%And here the minimum and maximum value of the colormap (0 for default):
caxis_min = 0;
caxis_max = 0;

%If you want to save the intensity map as png file, activate:
saveResult = 1;
%Notice: You have to create a folder 'Result' in this folder first!

%If you want to perform a 2D Gaussian fit of the intensity, activate:
performFit = 0;
%You can fit the whole image by using:
fitFullImage = 0;
%Or otherwise give the bounds for fitting in x and y here:
FIT_X_MIN = 102.2;
FIT_X_MAX = 103.2;
FIT_Y_MIN = 162.7;
FIT_Y_MAX = 163.7;
%If you want to see the residuals of this fit, activate:
plotResiduals = 0;

%Finally you have to set your step size and acquisition time of the scan:
PIXEL_SIZE       = 0.25;           %Effective pixel size in µm
ACQUISITION_TIME = 15;             %Acquisition time in seconds


%% ---- Data Readout ------------------------------------------------------

dataPath = ['Data\', dataFolder];
cd(dataPath)
data = importdata(fileName);        %Camera counts data
int = data(:,2:end);                %Cuts first column showing row number
int = int / ACQUISITION_TIME;       %Transforms intensity in cps
nPixels = length(int);              %Number of camera pixels
cd ..\..

%Define matrices for x- and y- position to create a mesh later:
x = zeros(nPixels);
y = zeros(nPixels);
for n = 1:nPixels
    x(:,n) = linspace(1,PIXEL_SIZE*nPixels,nPixels);
    y(n,:) = linspace(1,PIXEL_SIZE*nPixels,nPixels);
end

%Define maximum limits for the plot if full image is activated:
if generateFullImage
    PLOT_X_MIN = min(min(x)); %#ok<*UNRCH>
    PLOT_X_MAX = max(max(x));
    PLOT_Y_MIN = min(min(y));
    PLOT_Y_MAX = max(max(y));
end


%% ---- 2D Gaussian Fit ---------------------------------------------------

if performFit
    
    if fitFullImage
        xFit = x;
        yFit = y;
        intFit = int;
    else
        fitEdgesX = x(:,1)>FIT_X_MIN & x(:,1)<FIT_X_MAX;
        fitEdgesY = y(1,:)>FIT_Y_MIN & y(1,:)<FIT_Y_MAX;   
        xFit   = x  (fitEdgesX,fitEdgesY);
        yFit   = y  (fitEdgesX,fitEdgesY);
        intFit = int(fitEdgesX,fitEdgesY);
    end

    meanX = mean(xFit(:));
    meanY = mean(yFit(:));

    minInt = min(intFit(:));
    maxInt = max(intFit(:)); 
    
    %To perform the 2D Gauss fit, generate a mesh of positions:
    mesh = zeros(size(xFit,1),size(yFit,2),2);
    mesh(:,:,1) = xFit;
    mesh(:,:,2) = yFit;

    %The starting params and bounds for the fit must be limited by:
    guessParams = [maxInt,meanX,1,meanY,1,minInt];
    lb = [0.5*(maxInt-minInt),FIT_X_MIN,0,FIT_Y_MIN,0,minInt/2];
    ub = [2*(maxInt-minInt),FIT_X_MAX,10,FIT_Y_MAX,10,2*minInt];

    %Now perform the fit (fit function defined at the end of script):
    [fitResult,resnorm,residual,exitflag,output,lambda,jacobian]  = ...
        lsqcurvefit(@D2GaussFunction,guessParams,mesh,intFit,lb,ub);

    fitConfs = nlparci(fitResult,residual,'jacobian',jacobian,...
                                          'alpha',0.32);
    fitData  = D2GaussFunction(fitResult,mesh);

    %Extract the fit parameters:
    x_0           = fitResult(2);
    errX_0        = x_0-fitConfs(2,1);
    y_0           = fitResult(4);
    errY_0        = y_0-fitConfs(4,1);
    FWHM_x        = 2*sqrt(log(2)*2)*(fitResult(3));
    errFWHM_x     = FWHM_x-2*sqrt(log(2)*2)*(fitConfs(3,1));
    FWHM_y        = 2*sqrt(log(2)*2)*(fitResult(5));
    errFWHM_y     = FWHM_y-2*sqrt(log(2)*2)*(fitConfs(5,1));
    peakCounts    = fitResult(1);
    errPeakCounts = peakCounts-fitConfs(1,1);

end


%% ---- Generate intensity map --------------------------------------------

if generatePlots

    figure('Position', [200 100 1300 800])
    surf(x,y,int,'FaceAlpha',0.9)
    title('Camera Image')
    shading interp
    colormap(infernoColorMap)
    colorbar
    if caxis_min ~= caxis_max
        caxis([caxis_min,caxis_max])
    end
    xlabel('x-position in \mum')
    ylabel('y-position in \mum')
    zlabel('counts')
    if PLOT_DIM == 2
        view(2)
    end
    
    if performFit
        hold on
        surf(xFit,yFit,fitData,'LineStyle','-','FaceAlpha',0)
        hold off
    
        legend('off')
        dim = [0.25 0.3 0.4 0.5];
        str = {sprintf('stdev(x) = (%.0f +/- %.0f) nm',...
                        FWHM_x*1000/2.355, errFWHM_x*1000/2.355),...
           sprintf('stdev(y) = (%.0f +/- %.0f) nm',...
                        FWHM_y*1000/2.355, errFWHM_y*1000/2.355),...
           sprintf('cps(peak) = (%.2f +/- %.2f) kcps',...
                        peakCounts/1000,  errPeakCounts/1000)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
    
    set(gca,'FontSize',12)
    
end


%% ---- Save plot as png --------------------------------------------------

if saveResult
    cd 'Result'
    saveas(gcf, dataFolder, 'png');
    cd ..
end


%% ---- Plot Residuals ----------------------------------------------------

if plotResiduals
    
    figure
    surf(x,y,int-fitData,'FaceAlpha',0.9)
    xlim([PLOT_X_MIN,PLOT_X_MAX])
    ylim([PLOT_Y_MIN,PLOT_Y_MAX])
    title('Residuals')
    colormap(infernoColorMap)
    shading interp
    colorbar
    xlabel('x-position in \mum')
    ylabel('y-position in \mum')
    zlabel('counts')

end


%% ---- Used functions ----------------------------------------------------

function F = D2GaussFunction(x,xdata) %#ok<*DEFNU>
 F = x(1)*exp(-((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + ...
     (xdata(:,:,2)-x(4)).^2/(2*x(5)^2))) + x(6);
end


%% ---- Protocol of updates -----------------------------------------------

% 26.08.20 (V1.1):  - Selective analysis removed.
%
% 30.01.21 (V1.1):  - Warnings removed.
