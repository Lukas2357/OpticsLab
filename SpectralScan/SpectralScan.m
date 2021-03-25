%
%,------------------------------------------------------------------------,
%| Separate BG & Signal in spectral scan   (Lukas Hunold @01/02/21)  V1.2 |
%'------------------------------------------------------------------------'
%
%In this script several spectra of SiV centers are analysed. The background
%is estimated and in the spectral window of Raman and SiV peak the
%integrated bg and integrated signal are calculated.

close all
clear variables


%% ---- Optional functions ------------------------------------------------

%Name your asc files with a number from 0 to the last one and possibly add
%an extension of the name afterwards by activating:
fileExtension     = false;
fileExtensionName = 'whatever';
%Give the number of the last file name here:
nFiles = 80;
%If you want to plot the spectra, activate:
PlotSpectrum = 1;


%% ---- Parameters to be predefined ---------------------------------------

ACQUISITION_TIME = 1;       %acquisition time in seconds
STEPSIZE         = 1;       %stepsize of spectral scan in microns

RAMAN_WL_MIN = 715;         %minimum WL considered for Raman signal
RAMAN_WL_MAX = 725;         %maximum WL considered for Raman signal
SiV_WL_MIN   = 733.5;       %minimum WL considered for SiV signal
SiV_WL_MAX   = 746.5;       %maximum WL considered for SiV signal
Bg_WL_MIN    = 725;         %minimum WL considered for Bg averaging
Bg_WL_MAX    = 730;         %maximum WL considered for Bg averaging


%% ---- Data readout and calculation of signal and bg ---------------------

cd 'Data'
%Initialize vectors for signal and bg intesities:
IntSiV     = 1:nFiles;
IntSiVBg   = 1:nFiles;
IntRaman   = 1:nFiles;
IntRamanBg = 1:nFiles;
%Now go through all files and determine the bg by averaging the counts in
%the chosen bg window. Then sum up the correspoding counts in the SiV and 
%Raman window to get the bg in this regions. Afterwards subtract the bg
%from the total counts and calculate the signals by integrating the counts
%in the corresponding windows:
for file = 0:nFiles-1
    %Define current file:
    currentFile = sprintf('%d',file);
    %Load the data:
    if fileExtension
        Data = importdata([currentFile,fileExtensionName,'.asc']);
    else
        Data = importdata([currentFile,'.asc']);
    end
    %Print status of the program to console:
    fprintf('Loaded data file %.0f/%.0f \n',file+1,nFiles)
    %Set wavelength and signal arrays:
    WL  = Data(:,1);
    cps = Data(:,2)/ACQUISITION_TIME;
    %Plot the spectrum, if activated:
    if PlotSpectrum
        figure %#ok<*UNRCH>
        plot   (WL,cps,'k')
        title  (sprintf('Spectrum %.0f/%.0f',file+1,nFiles))
        xlabel ('WL / nm')
        ylabel ('cps')
        legend ('off')
        xlim   ([min(WL),max(WL)])
        ylim   ([0,max(cps)*1.1])
        set    (gca,'FOntSize',12)
        drawnow()
    end
    %Find the background baseline and subtract it from the data:
    bg  = mean(cps(WL<Bg_WL_MAX & WL>Bg_WL_MIN));
    cps = cps - bg;
    %And extract the different contributions:
    IntSiV(file+1)     = sum(cps(WL<SiV_WL_MAX & WL>SiV_WL_MIN));
    IntSiVBg(file+1)   = bg*length(WL(WL<SiV_WL_MAX & WL>SiV_WL_MIN));
    IntRaman(file+1)   = sum(cps(WL<RAMAN_WL_MAX & WL>RAMAN_WL_MIN));
    IntRamanBg(file+1) = bg*length(WL(WL<RAMAN_WL_MAX & WL>RAMAN_WL_MIN));
end
cd ..


%% ---- Plot Result -------------------------------------------------------

figure
posArray = 0:STEPSIZE:(nFiles-1)*STEPSIZE;
curve(1) = plot(posArray,IntRaman/1000);
hold on 
curve(2) = plot(posArray,IntSiV/1000);
curve(3) = plot(posArray,IntSiVBg/1000);
title  ('Spectral scan')
xlabel ('position on the sample / \mum')
ylabel ('kcps')
legend ('Raman intensity', 'SiV intensity', 'Bg around 738nm',...
        'Location', 'best')
xlim([0,(nFiles-1)*STEPSIZE])
set (gca,'FontSize',12)
set (curve,'linewidth',1)
    
    
%% ---- Protocol of updates -----------------------------------------------

% 26.08.20 (V1.1):  - Sections introduced.
%
% 01.02.21 (V1.2):  - Warnings removed.
