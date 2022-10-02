%
%,------------------------------------------------------------------------,
%| Calculating peak intensities of focal spots     (Lukas @27.01.21) V1.5 |
%'------------------------------------------------------------------------'
%
%In this script you can calculate the relative peak intensity of a focused 
%Gaussian beam for different filling factors (f0).

close all
clear variables


%% ---- Define parameters -------------------------------------------------

%First set here the parameters you use in the experiment:
NA = 0.95;            %Objective NA
P_INC = 7.8;          %Power before objective in mW
P_EFF = 5.1;          %Power after objective in mW
TRANS = 0.9;          %Fraction of transimtted power in the objective

%The corresponding objective filling is calculated out of that:
real_f0 = (2*log(P_INC/(P_INC-P_EFF/TRANS)))^(-1/2);

%For comparison also other fillings are simulated in the following, here
%you can choose which:
MIN_f0_SIM   = 0;       %minimum simulated filling factor
MAX_f0_SIM   = 6;       %minimum simulated filling factor
f0_STEPS_SIM = 100;     %number of steps simulated between the two


%% ---- Calculate relative peak intensity for different fillings ----------
    
f0    = real_f0;                   %First use the real objective filing
field = zeros(f0_STEPS_SIM,1);     %And initialize the field at the focus

%Now the fct for getting the focal field of a laser is used to get the
%laser intensity at the focus point for the given filling factor:
for f0Index=1:f0_STEPS_SIM
    %This is the fct under the integral, which gives the focus field:
    f   = @(x) (1-exp(-4./(2*f0.^2))).*NA^2*exp(-x.^2/f0.^2).*...   
                x.*((1-NA^2*x.^2).^(1/4)+(1-NA^2*x.^2).^(-1/4));
    %Integrate the fct to the max angle to get the field:
    field(f0Index) = integral(f,0,1);
    %Generate now new filling factors for next iteration:
    f0 = MIN_f0_SIM + (MAX_f0_SIM - MIN_f0_SIM)/f0_STEPS_SIM*f0Index;   
end

%Generate an array for the filling factors plottet in the following:
f0_array = linspace(MIN_f0_SIM,MAX_f0_SIM,f0_STEPS_SIM); 
%And calculate the intensity by squaring the field and normalizing it:
normIntensity = field.^2/max(field.^2); 
                                  
                                  
%% ---- Plot relative peak intensity for different fillings ---------------

figure 
plot(f0_array(2:end),normIntensity(2:end),'LineWidth',2)  
hold on
plot([min(f0_array),max(f0_array)],[normIntensity(1),normIntensity(1)],...
     'LineWidth',2)
title  ('Laser focus peak intensities depending on objective filling')
xlabel ('filling factor')
ylabel ('normalized intensity')
grid on 
legend off
str= ['P_{incident} = ',num2str(P_INC,2), 'mW \newline' , ...
      'P_{effective} = ',num2str(P_EFF,2), 'mW \newline \newline' , ...
      'Objective filling factor  = ',num2str(real_f0,2), '\newline', ...
      'Corresp. norm. Intens. = ',num2str(normIntensity(1),2)];
annotation('textbox',[0.5 0.76 0.07 0.02], 'String',{ str},...
           'FitBoxToText','on','BackgroundColor','w');

    
%% ---- Protocol of updates -----------------------------------------------

% 26.08.20 (V1.4):  - Sections introduced.
% 27.01.21 (V1.5):  - Warnings removed
%                   - Appearance and documentation improved
