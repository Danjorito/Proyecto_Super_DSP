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

denoised = zeros(length(wavelets), 4, 1000);

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

%%
%MOCAP Signal Replica

%Pre cleaning
close all;
clear;
clc;

%Loading the mocap signals from the original paper using the helper
%function amc_to_matrix.m

Mocap = amc_to_matrix("11_01.amc");

%Extraction of the elements used on the original paper
ch1 = Mocap(:, 2); %Called channel 1 in the original paper
ch1 = 6*(ch1-17.5)+86; %Replicating scaling and offset applied in the original paper
ch3 = Mocap(:, 6); %Called channel 3 in the original paper

%Figures of the extracted signals
figure(1)
plot(ch1);
title("Low frequency");
xlabel("Samples");
ylabel("Magnitude");

figure(2)
plot(ch3);
title("Low+High frequency");
xlabel("Samples");
ylabel("Magnitude");

%Uniform denoising using wavelets for ch1 and ch3

wavelets = ["haar", "db4", "db12", "sym2", "sym4", "sym8", "bl7", "bl10"];

denoised_ch1 = zeros(length(wavelets), 4, size(ch1, 1));
denoised_ch3 = zeros(length(wavelets), 4, size(ch3, 1));

for n = 1:size(denoised_ch1, 2) %This represents the levels of decomposition used on the denoising process
    for k = 1:size(denoised_ch1, 1)
        denoised_ch1(k, n, :) = wdenoise(ch1, n, Wavelet=wavelets(k));
        denoised_ch3(k, n, :) = wdenoise(ch3, n, Wavelet=wavelets(k));
    end
end

figure(3)
hold on;
plot(squeeze(denoised_ch3(3, 4, :)));
plot(ch3, 'LineWidth', 1.2);
legend("Sine wave with white Gaussian noise added", "Original signal");
hold off;

%ch1

%SM1
SM1_ch1 = zeros(size(denoised_ch1, 1), size(denoised_ch1, 2));
SM1_percentage_ch1 = zeros(size(denoised_ch1, 1), size(denoised_ch1, 2));

%In this case we calculate the SM1 for the base signal
signal_vector = ch1;
dif = signal_vector(2:end) - signal_vector(1:end-1);
SM1_base_ch1 = sum(abs(dif));

for n = 1:size(denoised_ch1, 2)
    for k = 1:size(denoised_ch1, 1)
        signal_vector = squeeze(denoised_ch1(k, n, :));
        dif = signal_vector(2:end) - signal_vector(1:end-1);
        SM1_ch1(k, n) = sum(abs(dif));
        SM1_percentage_ch1(k, n) = 100*(SM1_base_ch1-SM1_ch1(k, n))/(SM1_base_ch1);
    end
end

%SM2
SM2_ch1 = zeros(size(denoised_ch1, 1), size(denoised_ch1, 2));
SM2_percentage_ch1 = zeros(size(denoised_ch1, 1), size(denoised_ch1, 2));

signal_vector = ch1;
dif = signal_vector(1:end-2) - 2.*signal_vector(2:end-1) + signal_vector(3:end);
SM2_base_ch1 = sum((dif).^(2));

for n = 1:size(denoised_ch1, 2)
    for k = 1:size(denoised_ch1, 1)
        signal_vector = squeeze(denoised_ch1(k, n, :));
        dif = signal_vector(1:end-2) - 2.*signal_vector(2:end-1) + signal_vector(3:end);
        SM2_ch1(k, n) = sum((dif).^(2));
        SM2_percentage_ch1(k, n) = 100*(SM2_base_ch1-SM2_ch1(k, n))/(SM2_base_ch1);
    end
end

%ch3

%SM1
SM1_ch3 = zeros(size(denoised_ch3, 1), size(denoised_ch3, 2));
SM1_percentage_ch3 = zeros(size(denoised_ch3, 1), size(denoised_ch3, 2));

%In this case we calculate the SM1 for the base signal
signal_vector = ch3;
dif = signal_vector(2:end) - signal_vector(1:end-1);
SM1_base_ch3 = sum(abs(dif));

for n = 1:size(denoised_ch3, 2)
    for k = 1:size(denoised_ch3, 1)
        signal_vector = squeeze(denoised_ch3(k, n, :));
        dif = signal_vector(2:end) - signal_vector(1:end-1);
        SM1_ch3(k, n) = sum(abs(dif));
        SM1_percentage_ch3(k, n) = 100*(SM1_base_ch3-SM1_ch3(k, n))/(SM1_base_ch3);
    end
end

%SM2
SM2_ch3 = zeros(size(denoised_ch3, 1), size(denoised_ch3, 2));
SM2_percentage_ch3 = zeros(size(denoised_ch3, 1), size(denoised_ch3, 2));

signal_vector = ch3;
dif = signal_vector(1:end-2) - 2.*signal_vector(2:end-1) + signal_vector(3:end);
SM2_base_ch3 = sum((dif).^(2));

for n = 1:size(denoised_ch3, 2)
    for k = 1:size(denoised_ch3, 1)
        signal_vector = squeeze(denoised_ch3(k, n, :));
        dif = signal_vector(1:end-2) - 2.*signal_vector(2:end-1) + signal_vector(3:end);
        SM2_ch3(k, n) = sum((dif).^(2));
        SM2_percentage_ch3(k, n) = 100*(SM2_base_ch3-SM2_ch3(k, n))/(SM2_base_ch3);
    end
end