clc 
clear all
close all

t = 1/200:1/200:5;
fc = 2;

%Upper branch
It = [ones(1,200), ones(1,200), -ones(1,200), -ones(1,200), ones(1,200)];
ct = cos(2*pi*fc*t);
Itct = ct.*It;
figure
subplot(3,1,1);
plot(t,It);
subplot(3,1,2);
plot(t,ct);
subplot(3,1,3);
plot(t,Itct);

%Lower Branch
Qt = [-ones(1,200), ones(1,200), -ones(1,200), ones(1,200), ones(1,200)];
st = sin(2*pi*fc*t);
Qtst = Qt.*st;
figure
subplot(3,1,1);
plot(t,Qt);
subplot(3,1,2);
plot(t,st);
subplot(3,1,3);
plot(t,Qtst);

%QPSK modulated signal 
qpskst = Itct+(-Qtst);
figure
plot(t,qpskst);



