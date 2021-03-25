bin=res/1000;
clear LTa 
clear shift
clear t
clear tau1
 clear tau2
clear taumean
clear A1
clear A2
clear S
clear errTau
clear rerrTau
close all

Image2
cd(path2)
%addpath(path2);
options = optimoptions('lsqnonlin','Display', 'iter');
k=0; %index of y-direction
errTau=ones((size/step+1))*nan;
rerrTau=errTau;
tau1=ones((size/step+1))*nan;
while k<=size/step
    j=0; %index of x-direction
clear dat    
    
    
    
dat=importdata([num2str(k),'.txt']);


if mod(k,2)~=0            %corrects position of the histogram
    dat=dat(end:-1:1,:);
end

while j<=size/step
clear shift
clear t
clear fit

        j=j+1;
shift=dat(j,1:N);                             %reduction from the full histogram to the one which consists of real data. 

while shift(1)~=max(shift)           %negative delay is not posible in labview. Because of that a negative cylce is built in to the matlab code. 
shift=circshift(shift,1);
end
rshift=circshift(shift,floor(1.5/bin));
% shift=shift-rshift(1);
% rshift=rshift-rshift(1);

i=0;
 SNR=rms(shift(1:length(shift)))/(rms(rshift((length(rshift)-20):length(rshift))))
S(j,k+1)=SNR;
while shift(N-i) >1.10*min(shift)
        i=i+1;
end
shift(N-i:N)=min(shift);

%shift=shift/max(shift);
t=0:1:(length(shift)-1); 
t=t*bin;                                                                          %Time vector in ns 
%options.Algorithm = 'levenberg-marquardt';


if SNR>(9)&&SNR~=inf%10^(0.1*snr(rshift))>1
    xo=[(max(shift)-min(shift)),T1,mean(rshift((length(rshift)-20):length(rshift)))];

    n=length(shift);
    
[fit,resnorm,residual,exitflag,OUTPUT,LAMBDA,Jacobian]=lsqnonlin(@(p)(p(1).*exp(-1/p(2)*t)+p(3)-shift),xo,[0 0 0],[]);%,options);  %Lifetimefitting nonlinear least squares 
tau1(j,k+1)=fit(2);            %Lifetime from Fit
A1(j,k+1)= fit(1);             %Amplitude of the fit


%Estimating the error on the lifetime fit
 
Jacobian=full(Jacobian);
varp= resnorm*inv(Jacobian'*Jacobian)/n;
std=sqrt(diag(varp));
errTau(j,k+1)=std(2);
rerrTau(j,k+1)=std(2)/fit(2)*100;

% if fit(2)>0.9&&fit(2)<1.3%SNR>55
%        fit(2)
%     figure
%     plot(t,rshift)
%     hold on 
%     I=fit(1)*exp(-(t-1.5)/fit(2))+fit(3);
%     plot(t(floor(1.5/bin):length(t)),I(floor(1.5/bin):length(t)),'LineWidth',3)
%      hold off
%     xlabel('t in ns')
%     ylabel('counts')
%      name=['plot_x',num2str(j),'y',num2str(k)];
%     cd(path3)
%      saveas( gcf,name, 'png' );
%     cd(path2);
%     close(gcf)
%  
% end




%           if fit(1)<fit(3)                      %Filter if the backgroundsignal  is higher than the amplitude set lifetime of the point to zero
%               tau1(j,k+1)=0;
%           end
%   else 
%           A1(j,k+1)=0;
%           tau1(j,k+1)=0;
%           fit= [0,1,0];
end



%AT(j,k+1)=fit(1)*tau1(j,k+1);


LTa(j,k+1,:)=rshift(1:N);                 %writes the histogram in an 3d matrix for saving the acquired data sets for each pixel 
    
% if 10^(0.1*snr(rshift))<1
%               A1(j,k+1)=0;
%           tau1(j,k+1)=0;
%           fit(1)=0;
% end
j
k


    
end
    k=k+1
end

cd(path3)


save('Histogram.mat','LTa')
save('FLIM.mat','tau1')
figure('units','normalized','outerposition',[0 0 1 1])
imagesc(tau1)
colorbar 

view(0,-90)
saveas( gcf, 'lifetime', 'jpg' );
savefig('FLIM.fig')
figure('units','normalized','outerposition',[0 0 1 1])
surf(x,y,tau1)
colorbar 

view(90,-90)


title('Surf-plot')
shading interp
colorbar

xlabel('x-position in \mum')
ylabel('y-position in \mum')
zlabel('lifetime in ns')
sx=[min(min(x)), max(max(x))];
sy=[min(min(y)), max(max(y))];

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)

imagesc(sx,sy,intensity)
xlabel('y-position in \mum')
ylabel('x-position in \mum')
title('Intensity')
colorbar
view(0,-90)

subplot(2,2,2)

imagesc(sx,sy,tau1)
xlabel('y-position in \mum')
ylabel('x-position in \mum')
title('Lifetime')
h=colorbar;
ylabel(h,'\tau in ns')
view(0,-90)

subplot(2,2,3)

imagesc(sx,sy,S)
xlabel('y-position in \mum')
ylabel('x-position in \mum')
title('Signal to noise ratio')
h=colorbar;
ylabel(h,'SNR')
view(0,-90)


subplot(2,2,4)

imagesc(sx,sy,rerrTau)
xlabel('y-position in \mum')
ylabel('x-position in \mum')
title('relative error of the Lifetime')
h=colorbar;
xlabel(h,'realative error in %')
view(0,-90)


saveas( gcf, 'lifetimeimage', 'jpg' );

cd(path2)
%delete *.txt
