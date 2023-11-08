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

%MA WS=3, WT='sym'
MA_3_sym_signal = movmean(noisy_signal, 3);

%MA WS=5, WT='asym'
MA_5_sym_signal = movmean(noisy_signal, 5);

%MA WS=2, WT='asym'
MA_2_asym_signal = movmean(noisy_signal, [1 0]); %[1 0] means use only one element from the past and the current value for the avarage.

%MA WS=3, WT='asym'
MA_3_asym_signal = movmean(noisy_signal, [2 0]); %[1 0] means use two elements from the past and the current value for the avarage.

%Filtered noisy signal plots
figure(2)
subplot(2, 2, 1);
hold on;
plot(MA_2_asym_signal);
plot(original, 'LineWidth', 1.2);
legend("MA filtered signal", "Original signal");
title("a) WS=2, WT='asym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 2);
hold on;
plot(MA_3_asym_signal);
plot(original, 'LineWidth', 1.2);
legend("MA filtered signal", "Original signal");
title("a) WS=3, WT='asym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 3);
hold on;
plot(MA_3_sym_signal);
plot(original, 'LineWidth', 1.2);
legend("MA filtered signal", "Original signal");
title("a) WS=3, WT='sym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 4);
hold on;
plot(MA_5_sym_signal);
plot(original, 'LineWidth', 1.2);
legend("MA filtered signal", "Original signal");
title("a) WS=5, WT='sym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

%Performance metrics calculation

signal_vector = [noisy_signal' MA_2_asym_signal' MA_3_asym_signal' MA_3_sym_signal' MA_5_sym_signal']; %From this point onwards this order is used for the signals inside the vectors.

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