%
%,------------------------------------------------------------------------,
%| Compute refl. and transm. from Fresnel coeff.  (Lukas @30/01/21)  V1.3 |
%'------------------------------------------------------------------------'
%
%With this script you can calculate the reflection and transmission at an
%interface by using Fresnel coefficients. You just have to give the IORs.

close all
clear variables


%% ---- IORs --------------------------------------------------------------

n_air     = 1.000;
n_water   = 1.333;
n_ethanol = 1.360;
n_quartz  = 1.460;
n_pmma    = 1.490;
n_glass   = 1.520;
n_diamond = 2.409;


n1 = n_glass;
n2 = n_diamond;


%% ---- Calculate Reflection and Transmission -----------------------------

alphaDeg = 1:90;
alphaRad = alphaDeg/180*pi;

betaRad = asin(sin(alphaRad)*n1/n2);

r_s = - sin(alphaRad-betaRad)./sin(alphaRad+betaRad);
r_p = tan(alphaRad-betaRad)./tan(alphaRad+betaRad);

R_s = r_s.^2;
R_p = r_p.^2;

T_s = 1-R_s;
T_p = 1-R_p;


%% ---- Plot Reflection and Transmission ----------------------------------

figure('Position', [550 400 1000 350])
subplot(1,2,1)
plot(alphaDeg, R_s, 'LineWidth',1.5)
hold on
plot(alphaDeg, R_p, 'LineWidth',1.5)
xlabel ('angle of incidence in degrees')
ylabel ('Reflectivity')
legend ('R_s', 'R_p', 'Location', 'northwest')
set(gca,'FontSize',12)

subplot(1,2,2)
plot(alphaDeg, T_s, 'LineWidth',1.5)
hold on
plot(alphaDeg, T_p, 'LineWidth',1.5)
xlabel ('angle of incidence in degrees')
ylabel ('Transmitivity')
legend ('T_s', 'T_p', 'Location', 'northwest')
set(gca,'FontSize',12)


%% ---- Protocol of updates -----------------------------------------------

% 26.08.20 (V1.2):  - Sections introduced.
%
% 30.01.21 (V1.3):  - Subplots introduced.
