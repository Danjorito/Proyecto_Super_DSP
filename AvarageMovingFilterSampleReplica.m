%%
%Sine signal replica

%Pre cleaning
close all;
clear;
clc;

%Replicating sine signal
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
xlabel("Samples");
ylabel("Magnitude");
hold off;

%MA WS=3, WT='sym'
MA_3_sym_signal = movmean(noisy_signal, 3);

%MA WS=5, WT='sym'
MA_5_sym_signal = movmean(noisy_signal, 5);

%MA WS=2, WT='asym'
MA_2_asym_signal = movmean(noisy_signal, [1 0]); %[1 0] means use only one element from the past and the current value for the avarage.

%MA WS=3, WT='asym'
MA_3_asym_signal = movmean(noisy_signal, [2 0]); %[2 0] means use two elements from the past and the current value for the avarage.

%Filtered noisy signal plots
figure(2)
hold on;
plot(MA_2_asym_signal);
plot(original, 'LineWidth', 1.2);
legend("MA filtered signal", "Original signal");
title("a) WS=2, WT='asym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(3)
hold on;
plot(MA_3_asym_signal);
plot(original, 'LineWidth', 1.2);
legend("MA filtered signal", "Original signal");
title("b) WS=3, WT='asym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(4)
hold on;
plot(MA_3_sym_signal);
plot(original, 'LineWidth', 1.2);
legend("MA filtered signal", "Original signal");
title("c) WS=3, WT='sym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(5)
hold on;
plot(MA_5_sym_signal);
plot(original, 'LineWidth', 1.2);
legend("MA filtered signal", "Original signal");
title("d) WS=5, WT='sym'");
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


%ch1 MA signal filtering

%MA WS=3, WT='sym'
MA_3_sym_ch1 = movmean(ch1, 3);

%MA WS=5, WT='sym'
MA_5_sym_ch1 = movmean(ch1, 5);

%MA WS=2, WT='asym'
MA_2_asym_ch1 = movmean(ch1, [1 0]); %[1 0] means use only one element from the past and the current value for the avarage.

%MA WS=3, WT='asym'
MA_3_asym_ch1 = movmean(ch1, [2 0]); %[2 0] means use two elements from the past and the current value for the avarage.


%Plots for filtered signals on ch1

figure(3)
hold on;
plot(MA_2_asym_ch1, LineWidth=1.5);
plot(ch1, "k");
legend("MA filtered signal", "Original signal");
title("a) WS=2, WT='asym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(4)
hold on;
plot(MA_3_asym_ch1, LineWidth=1.5);
plot(ch1, "k");
legend("MA filtered signal", "Original signal");
title("b) WS=3, WT='asym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(5)
hold on;
plot(MA_3_sym_ch1, LineWidth=1.5);
plot(ch1, "k");
legend("MA filtered signal", "Original signal");
title("c) WS=3, WT='sym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(6)
hold on;
plot(MA_5_sym_ch1, LineWidth=1.5);
plot(ch1, "k");
legend("MA filtered signal", "Original signal");
title("d) WS=5, WT='sym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;


%Metrics for ch1 filtering

signal_vector = [ch1 MA_2_asym_ch1 MA_3_asym_ch1 MA_3_sym_ch1 MA_5_sym_ch1]; %From this point onwards this order is used for the signals inside the vectors.

%SM1
SM1_ch1 = zeros(1, size(signal_vector, 2));
SM1_percentage_ch1 = zeros(1, size(signal_vector, 2));
for n = 1:size(signal_vector, 2)
    dif_ch1 = signal_vector(2:size(signal_vector(:, n), 1), n) - signal_vector(1:size(signal_vector(:, n), 1)-1, n);
    SM1_ch1(n) = sum(abs(dif_ch1));
    SM1_percentage_ch1(n) = 100*(SM1_ch1(1)-SM1_ch1(n))/(SM1_ch1(1));
end

%SM2
SM2_ch1 = zeros(1, size(signal_vector, 2));
SM2_percentage_ch1 = zeros(1, size(signal_vector, 2));
for n = 1:size(signal_vector, 2)
    dif_ch1 = signal_vector(1:size(signal_vector(:, n), 1)-2, n) - 2.*signal_vector(2:size(signal_vector(:, n), 1)-1, n) + signal_vector(3:size(signal_vector(:, n), 1), n);
    SM2_ch1(n) = sum((dif_ch1).^(2));
    SM2_percentage_ch1(n) = 100*(SM2_ch1(1)-SM2_ch1(n))/(SM2_ch1(1));
end


%ch3 MA signal filtering

%MA WS=3, WT='sym'
MA_3_sym_ch3 = movmean(ch3, 3);

%MA WS=5, WT='sym'
MA_5_sym_ch3 = movmean(ch3, 5);

%MA WS=2, WT='asym'E_estE_est
MA_2_asym_ch3 = movmean(ch3, [1 0]); %[1 0] means use only one element from the past and the current value for the avarage.

%MA WS=3, WT='asym'
MA_3_asym_ch3 = movmean(ch3, [2 0]); %[2 0] means use two elements from the past and the current value for the avarage.


%Plots for filtered signals on ch3

figure(7)
hold on;
plot(MA_2_asym_ch3, LineWidth=1.5);
plot(ch3, "k");
legend("MA filtered signal", "Original signal");
title("a) WS=2, WT='asym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(8)
hold on;
plot(MA_3_asym_ch3, LineWidth=1.5);
plot(ch3, "k");
legend("MA filtered signal", "Original signal");
title("b) WS=3, WT='asym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(9)
hold on;
plot(MA_3_sym_ch3, LineWidth=1.5);
plot(ch3, "k");
legend("MA filtered signal", "Original signal");
title("c) WS=3, WT='sym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(10)
hold on;
plot(MA_5_sym_ch3, LineWidth=1.5);
plot(ch3, "k");
legend("MA filtered signal", "Original signal");
title("d) WS=5, WT='sym'");
xlabel("Samples");
ylabel("Magnitude");
hold off;


%Metrics for ch3 filtering

signal_vector = [ch3 MA_2_asym_ch3 MA_3_asym_ch3 MA_3_sym_ch3 MA_5_sym_ch3]; %From this point onwards this order is used for the signals inside the vectors.

%SM1
SM1_ch3 = zeros(1, size(signal_vector, 2));
SM1_percentage_ch3 = zeros(1, size(signal_vector, 2));
for n = 1:size(signal_vector, 2)
    dif_ch3 = signal_vector(2:size(signal_vector(:, n), 1), n) - signal_vector(1:size(signal_vector(:, n), 1)-1, n);
    SM1_ch3(n) = sum(abs(dif_ch3));
    SM1_percentage_ch3(n) = 100*(SM1_ch3(1)-SM1_ch3(n))/(SM1_ch3(1));
end

%SM2
SM2_ch3 = zeros(1, size(signal_vector, 2));
SM2_percentage_ch3 = zeros(1, size(signal_vector, 2));
for n = 1:size(signal_vector, 2)
    dif_ch3 = signal_vector(1:size(signal_vector(:, n), 1)-2, n) - 2.*signal_vector(2:size(signal_vector(:, n), 1)-1, n) + signal_vector(3:size(signal_vector(:, n), 1), n);
    SM2_ch3(n) = sum((dif_ch3).^(2));
    SM2_percentage_ch3(n) = 100*(SM2_ch3(1)-SM2_ch3(n))/(SM2_ch3(1));
end