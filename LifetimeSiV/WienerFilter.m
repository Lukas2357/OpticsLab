close all
clear variables

t = 0:1000;
tau1 = 100;
tau2 = 10;
signal = exp(-t/tau1);
irf = exp(-t/tau2);

figure
plot(t,signal)
hold on
plot(t,irf)

conv = ifft(fft(signal).*fft(irf));
figure
plot(t,conv)

conv = conv./max(conv);

%For using the Wiener-Filter, you have to determine the SNR.
    FFT_COEFFICIENT = 2^nextpow2(length(signal));
    signalFT = fft(conv,FFT_COEFFICIENT); 
    irfFT = fft(irf,FFT_COEFFICIENT);
    %Now the actual deconvolution is done by the filterfunction
    signalFT = (signalFT./irfFT); % .*...
        %((irfFT.*conj(irfFT))./((irfFT.*conj(irfFT))+1/SNR));
    %And the back FT
    signal = (ifft(signalFT));
    signal = signal(1:length(t));
    
figure
plot(t,signal)