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

%Uniform denoising using wavelets

wavelets = ["haar", "db4", "db12", "sym2", "sym4", "sym8", "bl7", "bl10"];

denoised = zeros(8, 4, 1000);

for n = 1:size(denoised, 2) %This represents the levels of decomposition used on the denoising process
    for k = 1:size(denoised, 1)
        denoised(k, n, :) = wdenoise(noisy_signal, n, Wavelet=wavelets(k));
    end
end

figure(2)

hold on;
plot(squeeze(denoised(1, 1, :)));
plot(original, 'LineWidth', 1.2);
legend("Sine wave with white Gaussian noise added", "Original signal");
hold off;

%SM1
SM1 = zeros(size(denoised, 1), size(denoised, 2));
SM1_percentage = zeros(size(denoised, 1), size(denoised, 2));

%In this case we calculate the SM1 for the base signal
signal_vector = noisy_signal';
dif = signal_vector(2:end) - signal_vector(1:end-1);
SM1_base = sum(abs(dif));

for n = 1:size(denoised, 2)
    for k = 1:size(denoised, 1)
        signal_vector = squeeze(denoised(k, n, :));
        dif = signal_vector(2:end) - signal_vector(1:end-1);
        SM1(k, n) = sum(abs(dif));
        SM1_percentage(k, n) = 100*(SM1_base-SM1(k, n))/(SM1_base);
    end
end

%SM2
SM2 = zeros(size(denoised, 1), size(denoised, 2));
SM2_percentage = zeros(size(denoised, 1), size(denoised, 2));

signal_vector = noisy_signal';
dif = signal_vector(1:end-2) - 2.*signal_vector(2:end-1) + signal_vector(3:end);
SM2_base = sum((dif).^(2));

for n = 1:size(denoised, 2)
    for k = 1:size(denoised, 1)
        signal_vector = squeeze(denoised(k, n, :));
        dif = signal_vector(1:end-2) - 2.*signal_vector(2:end-1) + signal_vector(3:end);
        SM2(k, n) = sum((dif).^(2));
        SM2_percentage(k, n) = 100*(SM2_base-SM2(k, n))/(SM2_base);
    end
end