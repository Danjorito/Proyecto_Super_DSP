%Pre cleaning
close all;
clear;
clc;

%Sine signal replica
t = 0:1:999;
original = sin(2*pi*0.002*t);

%Noisy replica
noise = wgn(1, 1000, -28); %-28 dBW was chosen to replicate the noise in the original paper.
noisy_signal = original + noise;
figure(1)
hold on;
plot(noisy_signal);
plot(original, 'LineWidth', 1.2);
legend("Sine wave with white Gaussian noise added", "Original signal");
hold off;

%Implementation of the Kalman filters, E_mea = 0.01
E_est = 0.01;
E_mea = 0.02;
EST = zeros(1, 1000);
EST(1) = noisy_signal(1);
KG = 0.0001;

for n = 1:size(noisy_signal, 2)-1
    EST(n+1) = EST(n);
    E_est = E_est + 0.01;
    EST(n+1) = EST(n+1) + KG*(noisy_signal(n+1) - EST(n+1));
    E_est = (1 - KG)*E_est;
    KG = (E_est)/(E_est + E_mea);
    disp("E_est: "+ num2str(E_est));
    disp("KG:    "+ num2str(KG));
end

figure(2)
hold on;
plot(EST);
plot(original, 'LineWidth', 1.2);
legend("Sine wave with white Gaussian noise added", "Original signal");
hold off;