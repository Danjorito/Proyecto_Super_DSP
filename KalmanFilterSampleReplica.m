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
E_est = 0.0016; %White noise standard deviation
E_mea = 0.01;
E_fut = 0.0003;
EST_1 = zeros(1, 1000);
EST_1(1) = noisy_signal(1); %Initial estimation

for n = 1:size(noisy_signal, 2)-1
    %Estimation part
    EST_1(n+1) = EST_1(n);
    E_est = E_est + E_fut;
    %Correction part
    KG = (E_est)/(E_est + E_mea);
    EST_1(n+1) = EST_1(n+1) + KG*(noisy_signal(n+1) - EST_1(n+1));
    E_est = (1 - KG)*E_est;
end

%Implementation of the Kalman filters, E_mea = 0.02
E_est = 0.0016;
E_mea = 0.02;
E_fut = 0.0001;
EST_2 = zeros(1, 1000);
EST_2(1) = noisy_signal(1);

for n = 1:size(noisy_signal, 2)-1
    %Estimation part
    EST_2(n+1) = EST_2(n);
    E_est = E_est + E_fut;
    %Correction part
    KG = (E_est)/(E_est + E_mea);
    EST_2(n+1) = EST_2(n+1) + KG*(noisy_signal(n+1) - EST_2(n+1));
    E_est = (1 - KG)*E_est;
    KG = (E_est)/(E_est + E_mea);
end

figure(2)
subplot(2, 1, 1)
hold on;
plot(EST_1);
plot(original, 'LineWidth', 1.2);
legend("Noisy sine wave with Kalman filter applied", "Original signal");
hold off;

subplot(2, 1, 2)
hold on;
plot(EST_2);
plot(original, 'LineWidth', 1.2);
legend("Noisy sine wave with  Kalman filter applied", "Original signal");
hold off;

%Performance metrics calculation

signal_vector = [noisy_signal' EST_1' EST_2']; %From this point onwards this order is used for the signals inside the vectors.

%SM1
SM1 = zeros(1, size(signal_vector, 2));
SM1_percentage = zeros(1, size(signal_vector, 2));
for n = 1:size(signal_vector, 2)
    dif = signal_vector(2:size(signal_vector(:, n), 1), n) - signal_vector(1:size(signal_vector(:, n), 1)-1, n);
    SM1(n) = sum(abs(dif));
    SM1_percentage(n) = 100*(SM1(1)-SM1(n))/(SM1(1));
end

%SM2
SM2 = zeros(1, size(signal_vector, 2));
SM2_percentage = zeros(1, size(signal_vector, 2));
for n = 1:size(signal_vector, 2)
    dif = signal_vector(1:size(signal_vector(:, n), 1)-2, n) - 2.*signal_vector(2:size(signal_vector(:, n), 1)-1, n) + signal_vector(3:size(signal_vector(:, n), 1), n);
    SM2(n) = sum((dif).^(2));
    SM2_percentage(n) = 100*(SM2(1)-SM2(n))/(SM2(1));
end