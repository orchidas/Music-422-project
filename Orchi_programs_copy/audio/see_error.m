close all;

[x,fs] = audioread('Castanets.wav');
[y,fs] = audioread('Castanets_bs_highDR.wav');

figure(1);
subplot(311);
plot(x)
subplot(312);
plot(y);
subplot(313);
plot(abs(x-y));