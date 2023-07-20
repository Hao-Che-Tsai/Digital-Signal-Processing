%Computing Assignment, Module 1.
%ECE 464/564 OSU

clc;
clear all;
close all;


%Generate two vectors with 1000 Gaussian random numbers each.
n_samples = 1000;
b1 = randn(1, n_samples);
b2 = randn(1, n_samples);

%Create two vectors of 1 and -1, where it is 1 if b1 >=0 and -1 otherwise.
b1(b1>=0) = 1;
b1(b1~=1) = -1;
b2(b2>=0) = 1;
b2(b2~=1) = -1;

%Upsampling with factor of 4
Up_factor1 = 4;
b1 = upsample(b1,Up_factor1);
b2 = upsample(b2,Up_factor1);

%implement the pulse shaping filter with 
Beta =  0.5;              %Roll off factor
sample_rate = 1;          %sampling rate
T = 4;                    
N = 32;
M=(-N/2):1:(N/2)-1;

% Sampling the frequency response of the square root raised cosine filter 
% at a high rate to avoid time-domain aliasing. The implementation assumes
% this length of H_rc is an even number.

H_rc = zeros(1,1000); 
i = 0;
for k= -length(H_rc)/2 : length(H_rc)/2-1  
    i = i+1;
    w = (2*pi*k)/length(H_rc);
    if abs(w) >= 0 && abs(w)<= pi*(1-Beta)/T
        H_rc(i) = sqrt(T);
    elseif abs(w) <= pi*(1+Beta)/T
        H_rc(i) = sqrt(0.5*T*(1+cos((T/(2*Beta))*(abs(w)-(pi*(1-Beta))/T))));
    end
    
end

% plotting the root raised cosine filter to make sure it is correct!
figure(1)
w = [-length(H_rc)/2:length(H_rc)/2-1]*2*pi/length(H_rc);
plot(w,H_rc)
ylabel('Amplitude','Fontsize',14,'FontWeight','Bold');
xlabel('Normalized Frequency (*2pi rad/s)','Fontsize',14,'FontWeight','Bold');
title('Frequency Response of Sqaure Root Raised Cosine Filter','Fontsize',8,'FontWeight','Bold');
set(gca,'Fontsize',11,'FontWeight','Bold','linewidth',.25); 
grid on

% Discrete Fourier transform requires H_rc to be defined in the range 
% [0,2*pi), Use the fftshift operation to accomplish this. 

H_rc = fftshift(H_rc);
w = fftshift(w); % Doing this to match the frequency correctly when plotting. 
w = 2*pi*(w<0)+ w; % Moving values of frequency in the range [-pi, 0) to the range [pi,2*pi)
figure(2)
plot(w,H_rc)
ylabel('Amplitude','Fontsize',14,'FontWeight','Bold');
xlabel('Normalized Frequency (*2pi rad/s)','Fontsize',14,'FontWeight','Bold');
title('Frequency Response of Sqaure Root Raised Cosine Filter','Fontsize',8,'FontWeight','Bold');
set(gca,'Fontsize',11,'FontWeight','Bold','linewidth',.25); 
grid on


%Take the inverse Fourier transform
h = ifft(H_rc);
%shift the second half of the inverse
h_rc = fftshift(h);

%Plot the time domain representation, h_rc[n]
figure(3)
n = -length(h_rc)/2:length(h_rc)/2-1;
stem(n,h_rc,'linewidth',1)
ylabel('Amplitude','Fontsize',14,'FontWeight','Bold');
xlabel('n','Fontsize',14,'FontWeight','Bold');
title(' Sqaure Root Raised Cosine Filter (Before Truncation)','Fontsize',8,'FontWeight','Bold');
set(gca,'Fontsize',14,'FontWeight','Bold','linewidth',.25); 
set(gcf, 'color', 'white');
grid on

% Pick the middle N samples to create an FIR filter with N coefficients

h_rc = h_rc(length(h_rc)/2 + 1 - N/2: length(h_rc)/2 +1 +N/2 - 1);

figure(4)
n = 0:N-1;
stem(n,h_rc,'linewidth',1)
ylabel('Amplitude','Fontsize',14,'FontWeight','Bold');
xlabel('n','Fontsize',14,'FontWeight','Bold');
title(' Sqaure root Raised Cosine Filter','Fontsize',8,'FontWeight','Bold');
set(gca,'Fontsize',14,'FontWeight','Bold','linewidth',.25); 
set(gcf, 'color', 'white');
grid on


%Apply the filter on generated sequences, b1 and b2.
x1 = fir_filter(b1,h_rc); 
x2 = fir_filter(b2,h_rc);

%Upsample the pulse shapped signals by a factor of 20

up_factor2 = 20;
x1_up = upsample(x1, up_factor2);
x2_up = upsample(x2, up_factor2);

%Implement the FIR low-pass filter
wc_cut = (3*pi/160);                 %lowpass filter cutoff freq
M = 500; %size of the window. This choice for the number of coefficients is arbitrary!
n = 0:M-1;

% Here we are calculating the inverse Fourier transform of the frequency
% response of the ideal lowpass filter. I have assumed that M is even and
% used a delay of M/2 rather than (M-1)/2 as was given in the handout. 
% (Both will work.)

for i = n
    if i-M/2 == 0 %(To avoid NaN)
        h_d(i+1) = 20*wc_cut/pi;
    else
        h_d(i+1) = 20*sin(wc_cut*(i-(M/2)))/(pi*(i-(M/2)));
    end
end

%Multiply the low-pass filter by a Hann window.
han = hann(length(h_d), 'symmetric');
h_d = h_d.*han';

%plot the frequency response of the filter.
figure(5);
subplot(2,1,1);
plot(h_d);
title('h_d: Impulse response function of the lowpass filter');
xlabel('time (samples)');ylabel('coefficient value');
H_n = fft(h_d);
H_n = fftshift(H_n);
subplot(2,1,2);
semilogy([-length(H_n)/2:length(H_n)/2-1]*2*pi/length(H_n),abs(H_n));
title('|H_d|: Magnitude response of the lowpass filter');xlabel('frequency (radians/sample)');ylabel('Magnitude');

%Pass the upsampled signal through the interpolation filter.
y1 = fir_filter(x1_up, h_d);
y2 = fir_filter(x2_up, h_d);


%Modulate the signals of the two paths:
w_c = 0.44*pi;                     %carrier frequency
t = 0:length(y1)-1;
mod1 = y1.*cos(w_c*t);
mod2 = y2.*sin(w_c*t);
%Merge the I and Q paths together
x_tran = mod1 + mod2;
n2=0:1:length(x_tran)-1;

figure(6)
plot(n2,x_tran)
ylabel('Amplitude','Fontsize',14,'FontWeight','Bold');
xlabel('Sample Number','Fontsize',14,'FontWeight','Bold');
title('Zoomed in Transmitted Signal','Fontsize',8,'FontWeight','Bold');
set(gca,'Fontsize',14,'FontWeight','Bold','linewidth',.25); 
set(gcf, 'color', 'white');
grid on
xlim([400 1000]);

% Plot the spectrum of clean signal; I use the matlab command pwelch here
% to estimate the spectrum. This is like taking the Fourier transform, but
% includes some averaging and gets cleaner results. 
figure(7) 
pwelch(x_tran)

% simulate the AWGN channel. 
% Generate a vector of zero-mean Gausssian random numbers. 
% Multiplying the noise samples by c makes the variance c^2. 
% Here, c = 0.1.
noise = 0.1*randn(1, length(x_tran));
%Add the AWGN noise to the pure signal.
s = x_tran + noise;
figure(7)
plot(n2,s)
ylabel('Amplitude','Fontsize',14,'FontWeight','Bold');
xlabel('Sample Number','Fontsize',14,'FontWeight','Bold');
title(' Zoomed in Transmitted Signal (with AWGN)','Fontsize',8,'FontWeight','Bold');
set(gca,'Fontsize',14,'FontWeight','Bold','linewidth',.25); 
set(gcf, 'color', 'white');
grid on
xlim([400 1000]);

% plot the spectrum of the noisy signal
figure(8)
pwelch(s)

function [y] = fir_filter(x,h)
h=h(:);
X=zeros(size(h));
y=zeros(size(x));
for n=1:length(x)
    X=[x(n);X(1:end-1)];
    y(n)=X'*h;
end
end

