
% Computer Project Module #2 
% close all
received_signal = s; % Comes from Module #1

upsample2 = 20; % upsample rate2 
downsample1 = 10; % A/D downsamplerate
downsample2 = 2;
omega_m = 1/upsample2*3/8; 

% Bandpass filtering% Design a FIR bandpass filter using a Hann window
order = 500; % Number of coefficients in BPF 
wc=0.44*pi; %carrier frequecy of received signal
wc_low = wc -omega_m*pi; % Lower cutoff frequency of BPF, H1
wc_high = wc + omega_m*pi; % Upper cutoff frequency of BPF, H1 

%generate the impulse response function of the bandpass filter
for i = 0: order-1
    if i-order/2 == 0
        h_BP(i+1) = (wc_high - wc_low)/pi;
    else
        h_BP(i+1) = sin(wc_high*(i-(order/2)))/(pi*(i-(order/2))) - sin(wc_low*(i-(order/2)))/(pi*(i-(order/2)));
    end
end

%Multiply the band-pass filter by a Hann window.
han = hann(length(h_BP), 'symmetric');
h_BP = h_BP.*han';

omi=-pi:(2/511*pi):pi; 
figure(),subplot(2,1,1),plot(omi,10*log10(abs(fft(h_BP,512)))),title('Frequency Response of BPF'); 
xlabel('omega');
ylabel('Gain in dB')
subplot(2,1,2),plot(omi,10*log10(abs(fft(received_signal,512))));
title('Frequency Spectrum of Received Signal'); 
xlabel('omega');
ylabel('Gain in dB')

r_BP = fir_filter(received_signal,h_BP);
figure(),
plot(r_BP);
title('Received Signal after BPF') 

%downsampling to perform A/D
r_downspl = downsample(r_BP,downsample1);
figure(),subplot(2,1,1),plot(r_downspl);
title('Signal afterdownsampling'),
subplot(2,1,2),plot(omi,10*log10(abs(fft(r_downspl,512)))),
title('Frequency Spectrum of Received Signal afterdownsampling');
xlabel('omega');
ylabel('Gain in dB') 

% Demodulation 
omega1 = 0.4*pi; % The carrier frequency after downsampling by 10 is 0.44*pi*10 = 4.4*pi, same as 0.4pi
r1_dem = r_downspl.* cos(omega1*[0:length(r_downspl)-1]); 
r2_dem = r_downspl.* sin(omega1*[0:length(r_downspl)-1]);
figure(), 
subplot(2,1,1), 
plot(r1_dem),
title('Demodulated Signal')
subplot(2,1,2),
stem(r2_dem); 
figure(),
subplot(2,1,1),
plot(omi,10*log10(abs(fft(r1_dem,512)))),
title('Frequency Spectrum of r1 after demodulation');
xlabel('omega');
ylabel('Gain in dB')
subplot(2,1,2),
plot(omi,10*log(abs(fft(r2_dem,512))));
title('FrequencySpectrum of r2 after demodulation');
xlabel('omega');
ylabel('Gain in dB') 

wc =  omega_m*downsample1*pi; % Cutoff frequency of the LPF, H2

order = 500; %size of the window. This choice for the number of coefficients is arbitrary!
for i = 0:order-1
    if i-order/2 == 0 %(To avoid NaN)
        h_d(i+1) = wc/pi;
    else
        h_d(i+1) = sin(wc*(i-(order/2)))/(pi*(i-(order/2)));
    end
end

%Multiply the low-pass filter by a Hann window.
han = hann(length(h_d), 'symmetric');
h_LP = h_d.*han';

r1_LP = fir_filter(r1_dem,h_LP);
r2_LP = fir_filter(r2_dem,h_LP);

% Downsampling by 2
r1_B = downsample(r1_LP,downsample2);
r2_B = downsample(r2_LP,downsample2);

figure(), 
subplot(2,1,1), 
plot(r1_B),
title('Signal at Point A and B, R1'); 
hold on; 
plot(x1,'r');
legend('Signal at Point B','Signal at Point A');
hold off 

subplot(2,1,2),
plot(r2_B),
title('Signal at Point A and B, R2');
hold on; 
plot(x2,'r');
legend('Signal at Point B','Signal at Point A');
hold off 

figure(), 
subplot(2,1,1), 
plot(r1_B),
title('Signal at Point A and B - Zoomed in, R1'); 
hold on; 
plot(x1,'r');
legend('Signal at Point B','Signal at Point A');
hold off 
xlim([0,500])

subplot(2,1,2),
plot(r2_B),
title('Signal at Point B - Zoomed in, R2');
hold on; 
plot(x2,'r');
legend('Signal at Point B','Signal at Point A');
hold off 
xlim([0,500])

function [y] = fir_filter(x,h)
h=h(:);
X=zeros(size(h));
y=zeros(size(x));
for n=1:length(x)
    X=[x(n);X(1:end-1)];
    y(n)=X'*h;
end
end
