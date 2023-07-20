%Computer HW 3
%Receiver Part B
%This project doesn't need a equlizer
close all;
%h_rc comes from module 1 and r1_B, r2_B come from module 2
%matched filter
h_mf=circshift(h_rc,31);
mf1=fir_filter(r1_B,h_mf);
mf2=fir_filter(r2_B,h_mf);

figure()
subplot(2,1,1);
plot(mf1);
title('b1 signal after matched filter');
subplot(2,1,2);
plot(mf2);
title('b2 signal after matched filter');

%h_d comes from module 1
%delay estimation
mf1=circshift(mf1,31);
mf2=circshift(mf2,31);
c_fd1=xcorr(h_d,mf1);
c_fd2=xcorr(h_d,mf2);

%downsample
c_fd1=downsample(c_fd1,4);
c_fd2=downsample(c_fd2,4);

figure()
subplot(2,1,1);
plot(c_fd1);
title('b1 signal after downsampling');
subplot(2,1,2);
plot(c_fd2);
title('b2 signal after downsampling');

%symbol dectaction
for nn=1:2000
    if c_fd1(nn)>=0
        b1n(nn)=1;
    elseif c_fd1(nn)<0
        b1n(nn)=-1;
    end
end
for nn=1:2000
    if c_fd2(nn)>=0
        b2n(nn)=1;
    elseif c_fd2(nn)<0
        b2n(nn)=-1;
    end
end

figure()
subplot(2,1,1);
plot(b1n);
title('b1 signal after symbol detection');
subplot(2,1,2);
plot(b2n);
title('b2 signal after symbol detection');


function [y] = fir_filter(x,h)
h=h(:);
X=zeros(size(h));
y=zeros(size(x));
for n=1:length(x)
    X=[x(n);X(1:end-1)];
    y(n)=X'*h;
end
end