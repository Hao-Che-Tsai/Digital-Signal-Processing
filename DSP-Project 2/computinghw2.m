clc;
close all;
clear all;

% parameters
bit_rate = 1000;            
sample = 8;    
fs = bit_rate*sample; %sampling frequency
fc1 = 2000; %carrier frequency
fc2 = 3000;
wc1=fc1/(fs/2);
wc2=fc2/(fs/2);
sourcen = 352; %signal length
rollof_f = 0.5; 

%random source
source = randi([0 1],1,sourcen);
fpre = ones(1,32); %synchronization
fbegin = [0 1 1 1 1 1 1 0]; %start frame
fhead = [fpre fbegin];  
fend = [0 1 1 1 1 1 1 0];   %end frame
fmsg = [fhead source fend];

%transmitter
bipolari_temp = zeros(1,length(fmsg));
bipolarq_temp = zeros(1,length(fmsg));
bipolari = fmsg(1:2:end);
bipolarq = fmsg(2:2:end);
for i = 1:length(bipolari)
    if(bipolari(i) ==1)
        bipolari_temp(2*i-1:2*i) = [1,1];
    else
        bipolari_temp(2*i-1:2*i) = [0,0];
    end
end
for i = 1:length(bipolarq)
    if(bipolarq(i) ==1)
        bipolarq_temp(2*i-1:2*i) = [1,1];
    else
        bipolarq_temp(2*i-1:2*i) = [0,0];
    end
end
bipolari = [bipolari_temp,0];
bipolarq = [0,bipolarq_temp];

%bipolar switching
bipolari = 2*bipolari-1;
bipolarq = 2*bipolarq-1;

%sampling of I
bipolarisource_temp = [bipolari',zeros(size(bipolari,2),sample-1)];
lengthx = size(bipolarisource_temp,1);
lengthy = size(bipolarisource_temp,2);
bipolarisource = reshape(bipolarisource_temp',1,lengthx*lengthy);

%sampling of Q
bipolarqsource_temp = [bipolarq',zeros(size(bipolarq,2),sample-1)];
lengthx = size(bipolarqsource_temp,1);
lengthy = size(bipolarqsource_temp,2);
bipolarqsource = reshape(bipolarqsource_temp',1,lengthx*lengthy);

%roll-off filter
rcosfir = rcosdesign(rollof_f,6,sample);
rcosisource = conv(bipolarisource,rcosfir);%channel I after sampling
rcosqsource = conv(bipolarqsource,rcosfir);%channel Q after sampling
filterdelay1 = (length(rcosfir)-1)/2;  %delay of roll-off filter

%modulate
%carrier trasmmiter
time = 1:length(rcosisource);
rcarrier = rcosisource.*cos(2*pi*fc1.*time/fs)-rcosqsource.*sin(2*pi*fc1.*time/fs);

%channel
ebn0 = -6:8;
snr = ebn0 - 10*log10(0.5*16);
t2=1:27072;
for i = 1:length(snr)
    %white gausse noise
    rcarrierwnoise=awgn(rcarrier,snr(i),'measured');

    %receiver demodulation 
    %bandpass filter
    firbp =fir1(128,[wc1 wc2],'bandpass',hamming(129));
    rcarrierwnoise=conv(firbp,rcarrierwnoise);

    %upsample
    rcarrierwnoise=upsample(rcarrierwnoise,8);

    %seperate signal with sine and cosine
    ri =rcarrierwnoise.*cos(2*pi*fc1.*t2/fs);
    rq =-(rcarrierwnoise.*sin(2*pi*fc1.*t2/fs));

    %lowpass filter
    firlp =fir1(128,0.4); 
    rlpi = conv(firlp,ri);
    rlpq = conv(firlp,rq);

    %downsample
    rsinidown = downsample(rlpi,2);
    rcosqdown = downsample(rlpq,2);
    %Figure 1 point B end at here
    
end

figure(1);
subplot(211);
plot(bipolari);
title('Time domain waveform of channel I');
subplot(212);
plot(bipolarq);
title('Time domain waveform of channel Q');

%waveform of channel I & Q after sampling
figure(2);
subplot(221);
plot(bipolarisource);
title('Time domain waveform of channel I');
subplot(222);
plot(abs(fft(bipolarisource)));
title('Frequency domain waveform of channel I');
subplot(223);
plot(bipolarqsource);
title('Time domain waveform of channel Q');
subplot(224);
plot(abs(fft(bipolarqsource)));
title('Frequency domain waveform of channel Q');

%waveform of channel I & Q after filter
figure(3);
subplot(221);
plot(rcosisource(1:1000));
title('Time domain channel I pass shaping filter');
subplot(222);
plot(abs(fft(rcosisource)));
title('Frequency domain channel I pass shaping filter');
subplot(223);
plot(rcosqsource(1:1000));
title('Time domain channel Q pass shaping filter');
subplot(224);
plot(abs(fft(rcosqsource)));
title('Frequency domain channel Q pass shaping filter');

%waveform after carrier modulation
figure(4);
subplot(211);
plot(rcarrier(1:500));
title('Time domain after carrier modulation');
subplot(212);
plot(abs(fft(rcarrier)));
title('Frequency domain after carrier modulation');

%signal after bandpass filter & upsampler
figure(5);
plot(rcarrierwnoise);
title('Modulation signal after bandpass filter & upsampler');

%Receiver signal in channel I & Q
figure(6);
subplot(211);
plot(ri);
title('Signal in channel I');
subplot(212);
plot(rq);
title('Signal in channel Q');

%Channel I & Q after lowpass filter
figure(7);
subplot(211);
plot(rlpi);
title('Channel I after lowpass filter');
subplot(212);
plot(rlpq);
title('Channel Q after lowpass filter');

%Channel I & Q after downsampler
figure(8);
subplot(211);
plot(rsinidown);
title('Channel I after downsampler');
subplot(212);
plot(rcosqdown);
title('Channel Q after downsampler');




