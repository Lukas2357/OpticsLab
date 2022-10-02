%
%,------------------------------------------------------------------------,
%| Calculate dipole collection efficiency (Lukas Hunold @17/01/21)  V 1.1 |
%'------------------------------------------------------------------------'
%
%In this script the collection efficiency of a dipole in a variable
%refrective index material using standard objectives is calculated. The
%dipole is assumed to be either horizontal, vertical or diagonal to the
%sample surface. Notice that the dipole needs to be located several
%wavelength below the surface, since otherwise near field effects occur.

close all 
clear variables


%% ---- Define initial parameters -----------------------------------------

n1 = 1;                      %IOR of the objective surrounding medium
n2 = 2.409;                  %IOR of the dipole surrounding medium
WL = 656*10^-9;              %Laser-Wavelength used for excitation
k1 = 2*pi/WL*n1;             %Absolute value of the k-vector in medium 1
k2 = 2*pi/WL*n2;             %Absolute value of the k-vector in medium 2
NA_range = 0.5:0.01:0.95;    %Range of numerical apertures simulated


%% ---- Calculate collection efficiencies ---------------------------------

%Initialize collection efficiency arrays for three dipole orientations:
nu_h = zeros(length(NA_range),1);
nu_v = zeros(length(NA_range),1);
nu_d = zeros(length(NA_range),1);

for iNA = 1:length(NA_range)

%Choose NA for present iteration and calculate max collection angle for it:
NA = NA_range(iNA);
theta_max = asin(NA/n2);

%Now calculate the collection efficiency by using the corresponding
%integrals over the emission pattern up to the maximum collection angle.

%Attention!!!: 
%The integration bounds and arguments of the sin functions 
%might not all be chosen correctly, review is needed here!

nu_h(iNA) = integral(@(theta) (2*k1*sqrt(1-sin(theta).^2)/...
                              (k1*sqrt(1-sin(theta).^2)*...
                              (n2/n1)+k2*sqrt(1-k1^2/k2^2*...
                              sin(theta).^2)*(n1/n2))).^2*sin(theta).^3*...
                              3/4,0,theta_max);
                        
nu_v(iNA) = integral(@(theta) (2*k1*sqrt(1-sin(theta+pi/2).^2)/...
                              (k1*sqrt(1-sin(theta+pi/2).^2)*...
                              (n2/n1)+k2*sqrt(1-k1^2/k2^2*...
                              sin(theta+pi/2).^2)*(n1/n2))).^2*...
                              sin(theta).^2.*sin(theta-theta_max+pi/2)*...
                              3/4,pi/2,pi/2+theta_max);
                           
nu_d(iNA) = integral(@(theta) (2*k1*sqrt(1-sin(theta).^2)/...
                              (k1*sqrt(1-sin(theta).^2)*...
                              (n2/n1)+k2*sqrt(1-k1^2/k2^2*...
                              sin(theta).^2)*(n1/n2))).^2*...
                              sin(theta+pi/4).^2.*sin(theta)*...
                              3/4,-theta_max,theta_max);

end


%% ---- Plot the results --------------------------------------------------

figure('Position',[1000,500,700,400])
plot(NA_range,nu_h*1000,'LineWidth',1.5)
hold on
plot(NA_range,nu_v*100,'LineWidth',1.5)
plot(NA_range,nu_d*100,'LineWidth',1.5)
legend('horizontal','vertical','diagonal','Location','NorthWest',...
       'FontSize',12)
xlabel('NA of collecting objective')
ylabel('Collection efficiency / %')
xlim([min(NA_range),max(NA_range)])
ylim([0,100*max([max(nu_v),max(nu_h),max(nu_d)])])
set(gca,'FontSize',12)


%% ---- Protocoll of updates ----------------------------------------------

% 17.01.21 (V1.1):  - Sections introduced
%                   - Fresnel coefficients squared
