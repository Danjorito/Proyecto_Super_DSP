%%
%Sine signal replica

%Pre cleaning
close all;
clear;
clc;

%Sine signal replication
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
hold on;
plot(EST_1);
plot(original, 'LineWidth', 1.2);
legend("Noisy sine wave with Kalman filter applied", "Original signal");
title("a) R = 0.01, Q = 0.0003");
hold off;

figure(3)
hold on;
plot(EST_2);
plot(original, 'LineWidth', 1.2);
legend("Noisy sine wave with  Kalman filter applied", "Original signal");
title("a) R = 0.02, Q = 0.0001");
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
ch1 = (6*(ch1-17.5)+86)'; %Replicating scaling and offset applied in the original paper
ch3 = Mocap(:, 6)'; %Called channel 3 in the original paper


% ch1 Kalman filtering

%Implementation of the Kalman filters, E_mea = 0.01
E_est = 0.0016; %White noise standard deviation
E_mea = 0.01;
E_fut = 0.00015;
EST_1_ch1 = zeros(1, size(ch1, 2));
EST_1_ch1(1) = ch1(1); %Initial estimation

for n = 1:size(ch1, 2)-1
    %Estimation part
    EST_1_ch1(n+1) = EST_1_ch1(n);
    E_est = E_est + E_fut;
    %Correction part
    KG = (E_est)/(E_est + E_mea);
    EST_1_ch1(n+1) = EST_1_ch1(n+1) + KG*(ch1(n+1) - EST_1_ch1(n+1));
    E_est = (1 - KG)*E_est;
end

%Implementation of the Kalman filters, E_mea = 0.02
E_est = 0.0016;
E_mea = 0.02;
E_fut = 0.00005;
EST_2_ch1 = zeros(1, size(ch1, 2));
EST_2_ch1(1) = ch1(1);

for n = 1:size(ch1, 2)-1
    %Estimation part
    EST_2_ch1(n+1) = EST_2_ch1(n);
    E_est = E_est + E_fut;
    %Correction part
    KG = (E_est)/(E_est + E_mea);
    EST_2_ch1(n+1) = EST_2_ch1(n+1) + KG*(ch1(n+1) - EST_2_ch1(n+1));
    E_est = (1 - KG)*E_est;
    KG = (E_est)/(E_est + E_mea);
end

%Implementation of the Kalman filters, E_mea = 0.005
E_est = 0.0016; %White noise standard deviation
E_mea = 0.005;
E_fut = 0.0002;
EST_3_ch1 = zeros(1, size(ch1, 2));
EST_3_ch1(1) = ch1(1); %Initial estimation

for n = 1:size(ch1, 2)-1
    %Estimation part
    EST_3_ch1(n+1) = EST_3_ch1(n);
    E_est = E_est + E_fut;
    %Correction part
    KG = (E_est)/(E_est + E_mea);
    EST_3_ch1(n+1) = EST_3_ch1(n+1) + KG*(ch1(n+1) - EST_3_ch1(n+1));
    E_est = (1 - KG)*E_est;
end

figure(1)
subplot(3, 1, 2)
hold on;
plot(EST_1_ch1, 'LineWidth', 1.2);
plot(ch1, "k");
legend("Noisy sine wave with Kalman filter applied", "ch1 signal");
title("b) Emea = 0.01");
hold off;

subplot(3, 1, 3)
hold on;
plot(EST_2_ch1, 'LineWidth', 1.2);
plot(ch1, "k");
legend("Noisy sine wave with  Kalman filter applied", "ch1 signal");
title("c) Emea = 0.02");
hold off;

subplot(3, 1, 1)
hold on;
plot(EST_3_ch1, 'LineWidth', 1.2);
plot(ch1, "k");
legend("Noisy sine wave with  Kalman filter applied", "ch1 signal");
title("a) Emea = 0.005");
hold off;


%Performance metrics calculation

signal_vector = [ch1' EST_3_ch1' EST_1_ch1' EST_2_ch1']; %From this point onwards this order is used for the signals inside the vectors.

%SM1
SM1_ch1 = zeros(1, size(signal_vector, 2));
SM1_percentage_ch1 = zeros(1, size(signal_vector, 2));
for n = 1:size(signal_vector, 2)
    dif = signal_vector(2:size(signal_vector(:, n), 1), n) - signal_vector(1:size(signal_vector(:, n), 1)-1, n);
    SM1_ch1(n) = sum(abs(dif));
    SM1_percentage_ch1(n) = 100*(SM1_ch1(1)-SM1_ch1(n))/(SM1_ch1(1));
end

%SM2
SM2_ch1 = zeros(1, size(signal_vector, 2));
SM2_percentage_ch1 = zeros(1, size(signal_vector, 2));
for n = 1:size(signal_vector, 2)
    dif = signal_vector(1:size(signal_vector(:, n), 1)-2, n) - 2.*signal_vector(2:size(signal_vector(:, n), 1)-1, n) + signal_vector(3:size(signal_vector(:, n), 1), n);
    SM2_ch1(n) = sum((dif).^(2));
    SM2_percentage_ch1(n) = 100*(SM2_ch1(1)-SM2_ch1(n))/(SM2_ch1(1));
end


% ch3 Kalman filtering

%Implementation of the Kalman filters, E_mea = 0.01
E_est = 0.0016; %White noise standard deviation
E_mea = 0.01;
E_fut = 0.00015;
EST_1_ch3 = zeros(1, size(ch3, 2));
EST_1_ch3(1) = ch3(1); %Initial estimation

for n = 1:size(ch3, 2)-1
    %Estimation part
    EST_1_ch3(n+1) = EST_1_ch3(n);
    E_est = E_est + E_fut;
    %Correction part
    KG = (E_est)/(E_est + E_mea);
    EST_1_ch3(n+1) = EST_1_ch3(n+1) + KG*(ch3(n+1) - EST_1_ch3(n+1));
    E_est = (1 - KG)*E_est;
end

%Implementation of the Kalman filters, E_mea = 0.02
E_est = 0.0016;
E_mea = 0.02;
E_fut = 0.00005;
EST_2_ch3 = zeros(1, size(ch3, 2));
EST_2_ch3(1) = ch3(1);

for n = 1:size(ch3, 2)-1
    %Estimation part
    EST_2_ch3(n+1) = EST_2_ch3(n);
    E_est = E_est + E_fut;
    %Correction part
    KG = (E_est)/(E_est + E_mea);
    EST_2_ch3(n+1) = EST_2_ch3(n+1) + KG*(ch3(n+1) - EST_2_ch3(n+1));
    E_est = (1 - KG)*E_est;
    KG = (E_est)/(E_est + E_mea);
end

%Implementation of the Kalman filters, E_mea = 0.005
E_est = 0.0016; %White noise standard deviation
E_mea = 0.005;
E_fut = 0.0002;
EST_3_ch3 = zeros(1, size(ch3, 2));
EST_3_ch3(1) = ch3(1); %Initial estimation

for n = 1:size(ch3, 2)-1
    %Estimation part
    EST_3_ch3(n+1) = EST_3_ch3(n);
    E_est = E_est + E_fut;
    %Correction part
    KG = (E_est)/(E_est + E_mea);
    EST_3_ch3(n+1) = EST_3_ch3(n+1) + KG*(ch3(n+1) - EST_3_ch3(n+1));
    E_est = (1 - KG)*E_est;
end

figure(2)
subplot(3, 1, 2)
hold on;
plot(EST_1_ch3, 'LineWidth', 1.2);
plot(ch3, "k");
legend("Noisy sine wave with Kalman filter applied", "ch3 signal");
title("b) Emea = 0.01");
hold off;

subplot(3, 1, 3)
hold on;
plot(EST_2_ch3, 'LineWidth', 1.2);
plot(ch3, "k");
legend("Noisy sine wave with  Kalman filter applied", "ch3 signal");
title("c) Emea = 0.02");
hold off;

subplot(3, 1, 1)
hold on;
plot(EST_3_ch3, 'LineWidth', 1.2);
plot(ch3, "k");
legend("Noisy sine wave with  Kalman filter applied", "ch3 signal");
title("a) Emea = 0.005");
hold off;


%Performance metrics calculation

signal_vector = [ch3' EST_3_ch3' EST_1_ch3' EST_2_ch3']; %From this point onwards this order is used for the signals inside the vectors.

%SM1
SM1_ch3 = zeros(1, size(signal_vector, 2));
SM1_percentage_ch3 = zeros(1, size(signal_vector, 2));
for n = 1:size(signal_vector, 2)
    dif = signal_vector(2:size(signal_vector(:, n), 1), n) - signal_vector(1:size(signal_vector(:, n), 1)-1, n);
    SM1_ch3(n) = sum(abs(dif));
    SM1_percentage_ch3(n) = 100*(SM1_ch3(1)-SM1_ch3(n))/(SM1_ch3(1));
end

%SM2
SM2_ch3 = zeros(1, size(signal_vector, 2));
SM2_percentage_ch3 = zeros(1, size(signal_vector, 2));
for n = 1:size(signal_vector, 2)
    dif = signal_vector(1:size(signal_vector(:, n), 1)-2, n) - 2.*signal_vector(2:size(signal_vector(:, n), 1)-1, n) + signal_vector(3:size(signal_vector(:, n), 1), n);
    SM2_ch3(n) = sum((dif).^(2));
    SM2_percentage_ch3(n) = 100*(SM2_ch3(1)-SM2_ch3(n))/(SM2_ch3(1));
end