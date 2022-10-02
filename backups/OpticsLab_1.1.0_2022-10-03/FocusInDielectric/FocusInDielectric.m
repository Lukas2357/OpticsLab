%
%,------------------------------------------------------------------------,
%| Focus a beam into a dielectric medium   (Lukas Hunold @30/01/21)  V2.6 |
%'------------------------------------------------------------------------'
%
%In this script the focusing of a laser beam into a dielectric medium is
%investigated. You might calculate the field strength along the beam axis
%and in radial direction for all three field components, as well as the
%corresponding intensity. You can set various incident parameters, like the
%beam profile, objective parameters, laser wavelength and IOR of the media.

%The situation is as follows: 
%An incident beam E_inc (default is a Gaussian beam profile, normalized 
%such that that the integral over the full x-y plane equals 1, with the 
%z-axis always referred to as the incident beam direction) is focused by an
%objective. 
%The beam is always linearly polarized along the x-axis (for other 
%polarizations you need to manualy modify the vector including the 
%Fresnel-coefficients). 
%The focusing is described by a transformation of the electric field vector
%into spherical coordinates and consideration of energy conservation. 
%The free parameters of the objective (NA, filling, working distance) will 
%determine the outcome of this focusing process. However you will not see 
%this outcome directly, because the next steps are applied immediately. 
%These are:
% 1) Refraction of the beam at the interface between the media, which is
%    described by applying Fresnel coefficients to the before transformed
%    electric field components
% 2) Calculation of the far field of the resulting beam in the second
%    medium with these coefficients
% 3) Calculation of the focused field by using the angular spectrum
%    representation
%All these steps can essentialy be sumarized into three single integrals
%(one for each field component) over the polar angle, up to the maximum 
%angle covered by the objective. (In practice the final substitution
%sin(theta) = q is done for simplification).
%You will find functions in the script called E_x, E_y, E_z, describing the
%functions to be integrated over. The arrays E_Tx, E_Ty and E_Tz are then
%the integrated fields, which are transmitted to the medium. You might
%calculate them for any rho, phi and z. In addition you can set z_0 for
%defining the position of the interface. The convention here is:
%The position of the focus without interface (n1=n2) is always at z=0. The
%distance of the interface from that is defined as z_0. Notice, that you
%should always put z_0<0 (otherwise you focus simply in n1), z_0>-WD
%(otherwise your interface is placed before the objective) and z>z_0
%(because you can only describe the transmitted field with the used
%integrals). Notice also, that your focus will be shifted to some position 
%z>0, if n2>n1 and to some z<0 in case n1>n2.
%
%   For more detailed information on all these calculations see
%   "Principles of Nano-Optics" (Novotny), p.56-66 and p.73-78

clear variables
close all


%% ---- Parameters to be predefined ---------------------------------------

n1 = 1;                  %IOR of medium 1 (surrounding the objective)
n2 = 2.409;              %IOR of medium 2 (where the laser is focused in)
NA = 0.95;               %NA of the objective
WD = 0.35*10^-3;         %Objective working distance
f  = sqrt(1-NA^2)*WD;    %Objective focal length     
f0 = 1;                  %Objective filling factor
WL = 656*10^-9;          %Laser-Wavelength
P  = 10^-3;              %Incident beam power in Watt
k1 = 2*pi/WL*n1;         %Absolute value of the k-vector in medium 1
k2 = 2*pi/WL*n2;         %Absolute value of the k-vector in medium 2


%% ---- Calculation and 3D plot of the focused field ----------------------
      
rho = linspace(0,5*WL,51);          %Simulated radii
phi = pi/4;                         %Simulated azimuthal angles
z   = linspace(-5*WL,5*WL,26);      %Simulated axial distances
z0  = linspace(0,-WL/2,2);          %Simulated interface positions
xyaxis = [fliplr(-rho)';rho']/WL;   %axis in x-y-plane in units of the WL
zaxis  = z/WL;                      %z-axis (both axes for visualization)

%Calculate maximum field and intensity (at focus, z0=0) for comparison:
E_max = squeeze(E_t(1,k1,k2,n1,n2,f0,NA,f,P,0,0,0,0));
I_max = abs(E_max).^2*10^-7;      %units: kW/cm^2

for i=1:length(z0)
 
    tic

    %Calculate the focal fields by the functions defined later:
    E_Tx = squeeze(E_t(1,k1,k2,n1,n2,f0,NA,f,P,z0(i),rho,phi,z));
    E_Ty = squeeze(E_t(2,k1,k2,n1,n2,f0,NA,f,P,z0(i),rho,phi,z));
    E_Tz = squeeze(E_t(3,k1,k2,n1,n2,f0,NA,f,P,z0(i),rho,phi,z));

    %Calculate the resulting intensities (for 'negative radii' also):
    Int_x = [flipud(abs(E_Tx).^2);abs(E_Tx).^2]*10^-7;
    Int_y = [flipud(abs(E_Ty).^2);abs(E_Ty).^2]*10^-7;
    Int_z = [flipud(abs(E_Tz).^2);abs(E_Tz).^2]*10^-7;

    %Plot intensity vs radius in 2D (figure outside plot function):
    nLines = 8;    %Number of lines shown in the plot (for different z)
    figure('Name',sprintf('Z_0 = %.1f WL',z0(i)/WL),...
           'Position', [0 500 900 500])
    Plot_2D(nLines,xyaxis,zaxis,I_max,Int_x,Int_y,Int_z)
    %And intensity vs radius and z in 3D (figure inside function):
    Plot_3D(0,xyaxis,zaxis,I_max,Int_x,Int_y,Int_z)

    fprintf('calculated radial focus %.0f/%.0f ',i,length(z0))
    
    toc
    pause(1)

end


%% ---- Focus position and FWHM depending on z_0 --------------------------
   
rho = 0;                            %Simulated radii
phi = pi/4;                         %Simulated azimuthal angles
z   = linspace(-10*WL,20*WL,2021);  %Simluated axial distances
z0  = linspace(0,-2*WL,51);         %Simulated interface positions
zaxis  = z/WL;                      %z-axis (both axes for visualization)
z0axis = z0/WL;                     %z0-axis (both axes for visualization)

zPeak = zeros(1,length(z0));        %initialize peak positions along z
FWHM  = zeros(1,length(z0));        %initialize peak FWHM along z

figure('Position', [0 100 900 300])

for i=1:length(z0)
    
    tic

    %Calculate field in x-direction (others 0 for rho=0):
    E  = squeeze(E_t(1,k1,k2,n1,n2,f0,NA,f,P,z0(i),rho,phi,z));
    %And the corresponding intensity:
    Int = abs(E).^2*10^-7;

    %Find peak position and FWHM of the intensity via fcts defined below:
    zPeak(i) = FindPeak(zaxis,Int)-z0(i)/WL;
    FWHM(i)  = FindFWHM(zaxis,Int);

    %Define number of lines shown in the intensity vs z-axis plot:
    nLines = 8;
    %Plot the intensity vs z-axis for some nLines equally spaced z0:
    zPlot(zaxis,Int,I_max,i,z0,nLines,WL)

    fprintf('calculated axial focus %.0f/%.0f ',i,length(z0))
    
    toc
    pause(0.1) 

end

%Plot and fit the FWHM depending on z0:
FWHMfitParams = fitFWHM(z0axis,FWHM);

%Plot and fit the peak position depending on z0:
PeakfitParams = fitPeak(z0axis,zPeak);


%% ---- Used functions ----------------------------------------------------

%The following function defines the transmission of the beam through the
%interface by applying Fresnel coefficient to the corresponding k-vector
%components for each direction. In addition you find Bessel- and sin/cos-
%functions, that are introduced when carrying out the Phi integration.

function fresnel = fresnel(dir,k1,k2,n1,n2,rho,phi,q)
    kz1 = k1*sqrt(1-q.^2);
    kz2 = k2*sqrt(1-k1^2/k2^2*q.^2);
    ts  = 2*kz1/(kz1+kz2);
    tp  = 2*kz1/(kz1*(n2/n1)+kz2*(n1/n2));
    x   = k1.*rho.*q;
    if dir == 1
        fresnel = besselj(0,x).*(ts+tp*kz2/k2)-...
              besselj(2,x).*(tp*kz2/k2-ts)*cos(2*phi);
    elseif dir == 2
        fresnel = - besselj(2,x).*(tp*kz2/k2-ts)*sin(2*phi);
    elseif dir == 3
        fresnel = 1i*besselj(1,x).*tp.*q/k2*k1*cos(phi);
    else
        fprintf('Use 1,2 or 3 as first argument for x, y or z direction')
    end
end

%By defining an incident field, multiplying it with above fresnel factors 
%and including phase shift in second medium, the infinity field E_inf can
%be calculated. This is then multiplied by additional terms in the so
%called angular spectrum representation (asr).

function asr = asr(dir,k1,k2,n1,n2,f0,NA,f,P,z0,rho,phi,z,q)
    kz1   = k1*sqrt(1-q.^2);
    kz2   = k2*sqrt(1-k1^2/k2^2*q.^2);
    %Incident field normalized to a total power of P:
    E_inc = sqrt(2*P/pi)/(f0*NA*f).*exp(-q.^2/(f0*NA)^2);
    E_inf = E_inc.*exp(1i.*(kz1-kz2)*z0).*...
            fresnel(dir,k1,k2,n1,n2,rho,phi,q);
    asr   = E_inf.*(1-q.^2).^(-1/4).*...
            k1/2.*q.*exp(1i.*kz2.*z)*1i*f*exp(-1i*k1*f);
end

%To get the transmitted field, just integrate the asr over the full polar
%angle covered by the objective. Here a substitution is made, such that
%sin(theta)=q and therefore the bounds are 0 and sin(theta_max)=NA.

function E_t = E_t(dir,k1,k2,n1,n2,f0,NA,f,P,z0,rho,phi,z)
E_t = zeros(length(rho),length(phi),length(z));
 for iRho = 1:length(rho)
  for iPhi = 1:length(phi)
   for iZ   = 1:length(z)
     E_t(iRho,iPhi,iZ) = integral(@(q)...
     asr(dir,k1,k2,n1,n2,f0,NA,f,P,z0,rho(iRho),phi(iPhi),z(iZ),q),0,NA);
   end
  end
 end
end

%To plot the intensity vs the xy-axis, the following function is used:

function Plot_2D(nLines,xyaxis,zaxis,I_max,Int_x,Int_y,Int_z)

    %Define array of indices for the z-values to be plottet:
    iZ = linspace(1,length(zaxis),nLines);
    
    for i=1:nLines
        
        %Round the index always towards z(i)=0, to get integers for the
        %indices to be plottet and have symmetric values around z=0:
        if zaxis(floor(iZ(i)))<0
        iZ(i) = ceil(iZ(i));
        else
        iZ(i) = floor(iZ(i));
        end
        
        subplot(2,2,1)
        Int = Int_x(:,iZ(i));
        plot(xyaxis,Int,'DisplayName',sprintf('z = %.1f WL',zaxis(iZ(i))));
        title('|E_x|^2')
        xlabel('radial distance / WL')
        ylabel('intensity / kW cm^{-2}')
        ylim([0,I_max])
        legend('-DynamicLegend');
        hold on
        
        subplot(2,2,2)
        Int = Int_y(:,iZ(i));
        plot(xyaxis,Int,'DisplayName',sprintf('z = %.1f WL',zaxis(iZ(i))));
        title('|E_y|^2')
        xlabel('radial distance / WL')
        ylabel('intensity / kW cm^{-2}')
        ylim([0,I_max/5000])
        legend('-DynamicLegend');
        hold on
        
        subplot(2,2,3)
        Int = Int_z(:,iZ(i));
        plot(xyaxis,Int,'DisplayName',sprintf('z = %.1f WL',zaxis(iZ(i))));
        title('|E_z|^2')
        xlabel('radial distance / WL')
        ylabel('intensity / kW cm^{-2}')
        ylim([0,I_max/300])
        legend('-DynamicLegend');
        hold on
        
        pause(1) 
    end
end

%To plot the intensity vs the xy- and z-axis in a surf-plot, the following
%function is used. You can provide the direction (1,2,3 for x,y,z) of the
%field, that you want to show as first argument. If you provide 0, all of
%them will be shown.

function Plot_3D(dir,xyaxis,zaxis,I_max,Int_x,Int_y,Int_z)

    p = inputParser;
    addOptional(p,'Int_x',0);
    addOptional(p,'Int_y',0);
    addOptional(p,'Int_z',0);
    parse(p,Int_x,Int_y,Int_z)

    mesh = zeros(length(xyaxis),length(zaxis),2);
    for i = 1:length(xyaxis)
        mesh(i,:,1) = zaxis;
    end
    for i = 1:length(zaxis)
        mesh(:,i,2) = xyaxis;
    end
    xyAxis = mesh(:,:,2);
    zAxis  = mesh(:,:,1);

    if dir ~= 0
        figure('Position', [950 500 960 500])
        if dir == 1
            Int = squeeze(p.Results.Int_x);
        elseif dir == 2
            Int = squeeze(p.Results.Int_y);
        elseif dir == 3
            Int = squeeze(p.Results.Int_z);
        end
            surf(xyAxis,zAxis,Int)
            colormap(infernoColorMap)
            xlabel('\rho / WL')
            ylabel('z / WL')
            zlim([0,I_max])
    else
        figure('Position', [950 500 960 500])
        
        subplot(2,2,1)
        Int = squeeze(p.Results.Int_x);
        surf(xyAxis,zAxis,Int)
        colormap(infernoColorMap)
        title('|E_x|^2')
        xlabel('\rho / WL')
        ylabel('z / WL')
        zlabel('intensity / kW cm^{-2}')
        zlim([0,I_max])
        
        subplot(2,2,2)
        Int = squeeze(p.Results.Int_y);
        surf(xyAxis,zAxis,Int*5000)
        title('|E_y|^2 x 5000')
        xlabel('\rho / WL')
        ylabel('z / WL')
        zlabel('intensity / kW cm^{-2}')
        zlim([0,I_max])
        
        subplot(2,2,3)
        Int = squeeze(p.Results.Int_z);
        surf(xyAxis,zAxis,Int*300)
        title('|E_z|^2 x 300')
        xlabel('\rho / WL')
        ylabel('z / WL')
        zlabel('intensity / kW cm^{-2}')
        zlim([0,I_max])
    end
end

%And to plot the focus along z for different z0, you can use:

function zPlot(zaxis,Int,I_max,i,z0,nLines,WL)
    %Define a z0 for each line, by picking equaly distant values from z0:
    iZ0 = floor(linspace(1,length(z0),nLines));
    %Plot the intensity vs z-axis curves for each picked z0:
    if ismember(i,iZ0)
        plot(zaxis,Int,'DisplayName',sprintf('z_0 = %.1f WL',z0(i)/WL));
        title('Focus along z for different z_0')
        xlabel('axial position / WL')
        ylabel('intensity / kW cm^{-2}')
        ylim([0,I_max])
        legend('-DynamicLegend');
        hold on
    end
end

%This function gives you the z-values corresponding to the peak position:
function Peak = FindPeak(axis,int)
    [~,imax] = max(int);
    Peak     = axis(imax);
end

%This function gives you the FWHM of the peaks along z:
function FWHM = FindFWHM(axis,int)
    FWHMaxis = axis(int>max(int)/2);
    FWHM  = (max(FWHMaxis)-min(FWHMaxis));
end

%This function plots and fits you the FWHM for different z0:
function fitPars = fitFWHM(z0axis,FWHM)
    f = fit(z0axis',FWHM','poly1');
    fitPars = [f.p1,f.p2];
    fitConfs  = confint(f);
    fitCurve = z0axis*f.p1+f.p2;
    figure('Position', [950 50 470 350])
    plot(z0axis,FWHM,'.k')
    hold on
    plot(z0axis,fitCurve,'-r')
    title('Axial FWHM vs z_0')
    xlabel('z_0 / WL')
    ylabel('FWHM / WL')
    ylim([0,max(FWHM)*1.5])
    annotation('textbox',[0.15 0.55 0.3 0.3],'String',...
    {sprintf('FWHM (z_0) = (%.2f +/- %.2f) z_0 + (%.2f +/- %.2f)',...
    fitPars(1),fitPars(1)-fitConfs(1,1),...
    fitPars(2),fitPars(2)-fitConfs(1,2))},'FitBoxToText','on');
end

%This function plots and fits you the peak position for different z0:
function fitPars = fitPeak(z0axis,zPeak)
    f = fit(z0axis',zPeak','poly1');
    fitPars = [f.p1,f.p2];
    fitConfs  = confint(f);
    fitCurve = z0axis*f.p1+f.p2;
    figure('Position', [1440 50 470 350])
    plot(z0axis,zPeak,'.k')
    hold on
    plot(z0axis,fitCurve,'-r')
    title('Axial focus position vs z_0')
    xlabel('z_0 / WL')
    ylabel('z_f / WL')
    ylim([0,max(zPeak)*1.5])
    annotation('textbox',[0.15 0.55 0.3 0.3],'String',...
    {sprintf('z_f (z_0) = (%.2f +/- %.2f) z_0 + (%.2f +/- %.2f)',...
    fitPars(1),fitPars(1)-fitConfs(1,1),...
    fitPars(2),fitPars(2)-fitConfs(1,2))},'FitBoxToText','on');
end


%% ---- Protocol of updates -----------------------------------------------

% 26.08.20 (V2.5):  - Sections introduced.
% 30.01.21 (V2.6):  - Warnings removed.
