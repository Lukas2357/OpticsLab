clear variables
close all

tau=2;                      %Estimates Lifetime (Stratingvalue)
bin=0.03;                   %Time resolution 0.004 for the picoharp300 with no binning
Tmax=50;                     %Maximal time given by the excitation rate

t=0:bin:Tmax;
t=t';
%D=importdata('Coumarin.txt');   %Example Data from PicoQuant for testing
%IRF=D(:,1);
%LTd=D(:,2);
IRF=importdata('Data\test_IRF.dat');   %IRF=instrument response function as a vector
%Enter here the name of your IRF without header. The time scale should be the same as for the measurement. 

LTd=importdata('Data\test.dat');    %LTd=lifetime histogram
%Enter here the name of your Lifetimedata.


IRF=IRF/max(IRF);
i=1;
while IRF(i)<1*max(IRF)      %Removes the first elemets of the IRF until the maximum. 
 i=i+1;
end
IRF(1:(i))=[];
IRF=IRF(1:length(t));
i=1;
while LTd(i)<1*max(LTd)      %Removes the first elemets of the IRF until the maximum. 
 i=i+1;
end
LTd(1:(i))=[];
LTd=LTd(1:length(t));
%LTd=LTd/max(LTd); 


%%%%%%%%%%%   Wiener-Filter%%%%%%%%%%%%%%%%%%%
    SNR=rms(LTd)/(rms(LTd((length(LTd)-200):length(LTd))))
      n=2^nextpow2(length(LTd));
      LTd=LTd-mean(LTd((length(LTd)-100):length(LTd)));
        G=fft(LTd,n); 
        
     H=fft(IRF(1:length(LTd)),n);
        %H(imag(H)>a)=a;
        %Filterfunction
        W=(1./H).*((H.*conj(H))./((H.*conj(H))+1/SNR));  
    F=(W.*G);
    
    %reconstructed signal
    f=(ifft(F));
    
    
    f=f(1:length(LTd));
    f(1)=[];
    
    %%%%%%% Fitting %%%%%%%%%%%%%%%%%%%%%%%%
    
    x0=[max(f)-min(f) tau min(f)];
    y0=[max(LTd)-min(LTd) tau min(LTd)];
    
fit2=lsqnonlin(@(d)(d(1).*exp(-d(2)^(-1).*(t))+d(3)-LTd),y0,[0 0 0]);   %Fit without deconvolution
t=t(1:(length(t)-1));
fit1=lsqnonlin(@(d)(d(1).*exp(-d(2)^(-1).*(t))+d(3)-f),x0,[0 0 0]);     %fit with deconvolution

disp('Lifetime without deconvolution')
T=fit2(2)
disp('Lifetime with deconvolution')
T=fit1(2)
