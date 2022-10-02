%
%,------------------------------------------------------------------------,
%| Plot and Fit confocal intensity scans   (Lukas Hunold @30/01/21)  V4.0 |
%'------------------------------------------------------------------------'
%
%In this script an Intensity map is generated from three txt files,
%including the x position (x.txt), y position (y.txt) and the intensity
%(int.txt). The files must be placed in a folder inside a folder 'Data' in
%this folder. Giving the name of this folder, the program will plot the
%data and perform a spatial Gaussian fit of the intensity, if needed.

close all
clear variables


%% ---- Optional functions ------------------------------------------------

%First give the folders inside folder 'Data', where the files are located:
dataFolders = {'test'};

%If you want to see both APDs, activate:
showBothAPDs      = 1;
%If you want to combine them (sum and difference), activate:
combineBothAPDs   = 1;
%If instead you want to see only one of it, give here the number (1 or 2):
showOnlyAPD       = 2;
%If you want to see the intensity map, activate:
generatePlots     = 1;
%You can plot them in one figure (they might overlap then), by activating
putScansTogether  = 1;
%If you want to shift the axes data, such that both start at 0, activate:
shiftAxisZero     = 1;
%If you want to see the full image in this case, activate:
generateFullImage = 1;
%Otherwise you have to give the x and y limits for the plot here:
PLOT_X_MIN = 87;   %#ok<*NASGU>
PLOT_X_MAX = 107;
PLOT_Y_MIN = 93;
PLOT_Y_MAX = 113;
%You can again shift the axes, if you still want to start at 0, by using:
shiftAxisZeroAgain = 0;
%Here you can set the plot dimension for visualization:
PLOT_DIM = 3;   %options: 2 for 2D plot, 3 for 3D plot
%You can normalize the intensity to the images mean value by activating:
normalizeCounts = 0;
%And here the minimum and maximum value of the colormap (0 for default):
caxis_min = 0;
caxis_max = 0;

%If you want to perform a 2D Gaussian fit of the intensity, activate:
performFit   = 1;
%You can fit the whole image by using:
fitFullImage = 1;
%Or otherwise give the bounds for fitting in x and y here:
FIT_X_MIN = 102.2;
FIT_X_MAX = 103.2;
FIT_Y_MIN = 162.7;
FIT_Y_MAX = 163.7;
%If you want to see the residuals of this fit, activate:
plotResiduals = 1;

%If you see a shift of subsequent lines in the scan, this might be due to a
%problem in data collection, which can be corrected by activating:
correctForLineShifting = 0;
%The program will shift the lines back by the following number of steps:
LINE_SHIFT = 1;

%You can bin the image by activating:
activateBinning = 0;
%And give the number of pixels to be binned in row and column:
L_BIN = 1;

%Finally you have to set your step size and acquisition time of the scan:
STEP_SIZE        = 0.1;           %Stepsize of the scan in micrometer
ACQUISITION_TIME = 1;             %Acquisition time in seconds


%% ---- Data Readout and Ordering -----------------------------------------

%Initialize the cell arrays for data storage:
x      = cell(length(dataFolders),1);        %Initialize x-Position data  
y      = cell(length(dataFolders),1);        %Initialize y-Position data
apd    = cell(length(dataFolders),2);        %Initialize both APD counts
binCps = cell(length(dataFolders),2);        %Initialize binned APD cps
xBin   = cell(length(dataFolders),1);        %Initialize binned x
yBin   = cell(length(dataFolders),1);        %Initialize binned y
maxCps = zeros(length(dataFolders),2);       %Maximum counts on each APD

%And start data reading and ordering:
for file = 1:length(dataFolders)
    
    %Move to data folder and load the matrices from the txt files:
    dataPath = ['Data\', dataFolders{file}];
    cd(dataPath)
    x{file}     = importdata('x.txt');       %Load x-Position data     
    y{file}     = importdata('y.txt');       %Load y-Position data
    apd{file,1} = importdata('int0.txt');    %Load counts on APD 1
    apd{file,2} = importdata('int.txt');     %Load counts on APD 2
    cd ..\..
    
    %If activated shift the position axes to 0:
    if shiftAxisZero
        if file == 1
            minX = min(x{1}(:,1));
            minY = min(y{1}(1,:));
        end
        x{file} = x{file} - minX;         %x-Position data     
        y{file} = y{file} - minY;         %y-Position data 
    end
    
    %The scans can show shifting of subsequent lines, caused by the piezo.
    %This is a bug in data collection and can be fixed by shifting each 2nd
    %line by a few steps. This is done in the following:
    if correctForLineShifting
        for iX = 1:length(x{file}(1,:)) %#ok<*UNRCH>
            if mod(iX,2) == 0
                x{file}(:,iX) = x{file}(:,iX) + LINE_SHIFT*STEP_SIZE;
            end
        end
    end
    
    for iAPD = 1:2
        
        %If activated normalize the image counts to the maximum:
        if normalizeCounts
            apd{file,iAPD} = apd{file,iAPD}/mean(mean(apd{file,iAPD}));
        end
        
        %Transform counts in cps by deviding through acq. time:
        apd{file,iAPD} = apd{file,iAPD}/ACQUISITION_TIME;
        %And record the maximum cps for normalization later:
        maxCps(file,iAPD) = max(max(apd{file,iAPD}));
        
        %If activated, the picture is binned by the following code:
        if activateBinning
            nBINSx = floor(length(x{file}(:,1))/L_BIN);
            nBINSy = floor(length(y{file}(1,:))/L_BIN);
            binCps{file,iAPD} = zeros(nBINSx,nBINSy);
            xBin{file}        = zeros(nBINSx,nBINSy);
            yBin{file}        = zeros(nBINSx,nBINSy);
            for iX=1:nBINSx
                for iY=1:nBINSy
                    bin = apd{file,iAPD}((iX-1)*L_BIN+1:iX*L_BIN,...
                                         (iY-1)*L_BIN+1:iY*L_BIN);              
                    binCps{file,iAPD}(iX,iY) = mean(bin(:));
                    xBin{file}(iX,iY) = x{file}(L_BIN*iX,L_BIN*iY);
                    yBin{file}(iX,iY) = y{file}(L_BIN*iX,L_BIN*iY);
                end
            end
            apd{file,iAPD} = binCps{file,iAPD};
            x{file} = xBin{file};
            y{file} = yBin{file};
        end
    end
    
    %If activated, combine the two APDs as difference and sum:
    if combineBothAPDs
        apd{file,1} = (apd{file,1} - apd{file,2});
        apd{file,2} = (apd{file,1} + apd{file,2});
    end
    
end


%% ---- 2D Gaussian Fit ---------------------------------------------------

if performFit
    
    %Initialize cell arrays for fit params:
    x0         = cell(length(dataFolders),2);
    x0std      = cell(length(dataFolders),2);
    y0         = cell(length(dataFolders),2);
    y0std      = cell(length(dataFolders),2);
    xFWHM      = cell(length(dataFolders),2);
    xFWHMstd   = cell(length(dataFolders),2);
    yFWHM      = cell(length(dataFolders),2);
    yFWHMstd   = cell(length(dataFolders),2);
    peakCps    = cell(length(dataFolders),2);
    peakCpsStd = cell(length(dataFolders),2);
    fitData    = cell(length(dataFolders),2);
    xFit       = cell(length(dataFolders),1);
    yFit       = cell(length(dataFolders),1);
    apdFit     = cell(length(dataFolders),2);
    
    %Perform fit for each file and APD separately:
    for file = 1:length(dataFolders)
        for iAPD = 1:2
            %Set the fit bounds in x and y:
            if fitFullImage
                xFit{file}        = x{file};
                yFit{file}        = y{file};
                apdFit{file,iAPD} = apd{file,iAPD};
            else
                fitEdgesX = y{file}(1,:)>FIT_X_MIN & ...
                            y{file}(1,:)<FIT_X_MAX;
                fitEdgesY = x{file}(:,1)>FIT_Y_MIN & ...
                            x{file}(:,1)<FIT_Y_MAX;   
                xFit{file}        = y{file}  (fitEdgesY,fitEdgesX);
                yFit{file}        = x{file}  (fitEdgesY,fitEdgesX);
                apdFit{file,iAPD} = apd{file,iAPD}(fitEdgesY,fitEdgesX);
            end
            
            %Extract values for fit starting params:
            minX   = min(xFit{file}(:));
            maxX   = max(xFit{file}(:));
            minY   = min(yFit{file}(:));
            maxY   = max(yFit{file}(:));
            meanX  = mean(xFit{file}(:));
            meanY  = mean(yFit{file}(:));
            minCps = min(apdFit{file,iAPD}(:));
            maxCps = max(apdFit{file,iAPD}(:)); 

            %To perform the 2D Gauss fit, generate a mesh of positions:
            mesh = zeros(size(xFit{file},1),size(yFit{file},2),2);
            mesh(:,:,1) = xFit{file};
            mesh(:,:,2) = yFit{file};

            %The starting params and bounds for the fit are limited by:
            guessParams = [maxCps,meanX,1,meanY,1,minCps];
            lb = [0.5*(maxCps-minCps),minX,0,minY,0,minCps];
            ub = [2.5*(maxCps-minCps),maxX,10*(maxX-minX), ...
                                      maxY,10*(maxY-minY), maxCps];

            %Now perform the fit (fit function defined below):
            [fitRes,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                    lsqcurvefit(@D2GaussFunction,guessParams,...
                                mesh,apdFit{file,iAPD},lb,ub);

            %And extract the parameters as well as the fit data:
            fitConfs = nlparci(fitRes,residual,'jacobian',jacobian,...
                                               'alpha',0.32);
            fitData{file,iAPD}  = D2GaussFunction(fitRes,mesh);

            %Extract the concrete fit parameters with uncertainties:
            x0{file,iAPD}         = fitRes(2);
            x0std{file,iAPD}      = x0{file,iAPD} - fitConfs(2,1);
            y0{file,iAPD}         = fitRes(4);
            y0std{file,iAPD}      = y0{file,iAPD} - fitConfs(4,1);
            xFWHM{file,iAPD}      = 2*sqrt(log(2)*2)*(fitRes(3));
            xFWHMstd{file,iAPD}   = xFWHM{file,iAPD} - 2*sqrt(log(2)*2)*...
                                    (fitConfs(3,1));
            yFWHM{file,iAPD}      = 2*sqrt(log(2)*2)*(fitRes(5));
            yFWHMstd{file,iAPD}   = yFWHM{file,iAPD} - 2*sqrt(log(2)*2)*...
                                    (fitConfs(5,1));
            peakCps{file,iAPD}    = fitRes(1);
            peakCpsStd{file,iAPD} = peakCps{file,iAPD} - fitConfs(1,1);
        end
    end
end


%% ---- Generate intensity map --------------------------------------------

if generatePlots
    
    %If the image bouds in x and y are given above, all plots will be
    %generated with those, except the case that generateFullImage is
    %activated. Then all will be shown fully. If putScansTogether is
    %activated, the same is true, but the all scans are put in one plot:
    if generateFullImage
        PLOT_X_MIN = min(x{1}(:));
        PLOT_X_MAX = max(x{1}(:));
        PLOT_Y_MIN = min(y{1}(:));
        PLOT_Y_MAX = max(y{1}(:));
        for file = 2:length(dataFolders)
            if min(x{file}(:)) < PLOT_X_MIN
                PLOT_X_MIN = min(x{file}(:));
            end
            if max(x{file}(:)) > PLOT_X_MAX
                PLOT_X_MAX = max(x{file}(:));
            end
            if min(y{file}(:)) < PLOT_Y_MIN
                PLOT_Y_MIN = min(y{file}(:));
            end
            if max(y{file}(:)) > PLOT_Y_MAX
                PLOT_Y_MAX = max(y{file}(:));
            end
        end
    end
    
    %If the scans are put together create just one figure window:
    if putScansTogether
        figure('Position', [150 150 1450 750])
    end
    %Then iterate over all files:
    for file = 1:length(dataFolders)
        %If the scans are not put together, create a window for each:
        if ~putScansTogether
                figure('Position', [150 150 1450 750])
        end
        %Then iterate over the two APDs:
        for iAPD = 1:2
            %Create a subplot for each APD selected to be shown:
            if showBothAPDs
                if performFit && plotResiduals 
                    subplot(2,2,iAPD) 
                else
                    subplot(1,2,iAPD)
                end
            else
                if performFit && plotResiduals 
                    subplot(2,1,1) 
                else
                    subplot(1,1,1)
                end
                if showOnlyAPD == 1
                    iAPD = 1; %#ok<FXSET>
                else
                    iAPD = 2; %#ok<FXSET>
                end
            end
            %And fill the scans together with the fits (if activated):
            surf(x{file},y{file},apd{file,iAPD},'FaceAlpha',0.9)
            if PLOT_DIM == 2
                view(2)
            end
            xlim([PLOT_X_MIN,PLOT_X_MAX])
            ylim([PLOT_Y_MIN,PLOT_Y_MAX])
            if caxis_min ~= caxis_max
                caxis([caxis_min,caxis_max])
            end
            if ~combineBothAPDs
                title(sprintf('Confocal scan APD %.0f',iAPD))
            else
                if iAPD == 1
                    title(sprintf('Difference of the APDs'))
                else
                    title(sprintf('Sum of the APDs'))
                end
            end
            shading interp
            colormap(infernoColorMap)
            hold on
            if performFit
                surf(xFit{file},yFit{file},fitData{file,iAPD}, ...
                    'LineStyle','-','FaceAlpha',0)
                legend('off')
                dim = [0.35+(iAPD-1)*0.45 0.4 0.4 0.5];
                str = {sprintf('std(x) = (%.0f +/- %.0f) nm',...
                                xFWHM{file,iAPD}*1000/2.355, ...
                                xFWHMstd{file,iAPD}*1000/2.355),...
                       sprintf('std(y) = (%.0f +/- %.0f) nm',...
                                yFWHM{file,iAPD}*1000/2.355, ...
                                yFWHMstd{file,iAPD}*1000/2.355),...
                       sprintf('cps    = (%.2f +/- %.2f) kcps',...
                                peakCps{file,iAPD}/1000,  ...
                                peakCpsStd{file,iAPD}/1000)};
                annotation('textbox',dim,'String',str,'FitBoxToText','on');
                %If activated make subplots for residuals and plot them:
                if plotResiduals
                    if showBothAPDs
                        subplot(2,2,iAPD+2)
                    else
                        subplot(2,1,2)
                    end
                    surf(xFit{file},yFit{file},...
                         fitData{file,iAPD}-apd{file,iAPD},'FaceAlpha',0.9)
                    colormap(infernoColorMap)
                    xlim([PLOT_X_MIN,PLOT_X_MAX])
                    ylim([PLOT_Y_MIN,PLOT_Y_MAX])
                    title(sprintf('Residuals of fit %.0f',iAPD))
                    shading interp
                    colorbar
                    xlabel('x-position in \mum')
                    ylabel('y-position in \mum')
                    zlabel('(fit - data) / cps')
                end
            end
            shading interp
            colormap(infernoColorMap)
            colorbar
            xlabel('x-position in \mum')
            ylabel('y-position in \mum')
            xlim([PLOT_X_MIN,PLOT_X_MAX])
            ylim([PLOT_Y_MIN,PLOT_Y_MAX])
            zlabel('cps')
            if PLOT_DIM == 2
                view(2)
            end
        end
    end
    
end


%% ---- Used functions ----------------------------------------------------

function F = D2GaussFunction(x,xdata) %#ok<*DEFNU>
 F = x(1)*exp(-((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + ...
     (xdata(:,:,2)-x(4)).^2/(2*x(5)^2))) + x(6);
end


%% ---- Protocol of updates -----------------------------------------------

% 05.02.20 (V3.6):  - cd 'Result' inside saveResult if statement
%                   - Introduction of sections for different parts
%                   - Circle function removed (was unused)
%                   - caxis limit option introduced
%
% 07.02.20 (V3.7):  - Selective analysis option added
%
% 26.08.20 (V3.8):  - Option for APD averaging added
%                   - Option for normalizing counts added
%                   - Option for multiple input files added
%
% 02.09.20 (V3.9):  - Problem of fit bounds solved
%
% 02.09.20 (V4.0):  - Inferno color map added
%                   - cell arrays introduced
%                   - putScansTogether option introduced
%                   - Selective analysis option removed
%                   - Appearance and documentation improved
