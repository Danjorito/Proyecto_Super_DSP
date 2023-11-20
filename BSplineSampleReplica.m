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

%Spline function fitting using the spap2 function
spline_app_500ctrl = spap2(500, 3, t, noisy_signal);
spline_app_200ctrl = spap2(200, 3, t, noisy_signal);
spline_app_100ctrl = spap2(100, 3, t, noisy_signal);
spline_app_50ctrl = spap2(50, 3, t, noisy_signal);

fn_form_500ctrl = fnval(fn2fm(spline_app_500ctrl, 'ppform'), t);
fn_form_200ctrl = fnval(fn2fm(spline_app_200ctrl, 'ppform'), t);
fn_form_100ctrl = fnval(fn2fm(spline_app_100ctrl, 'ppform'), t);
fn_form_50ctrl = fnval(fn2fm(spline_app_50ctrl, 'ppform'), t);

figure(2)
hold on;
plot(fn_form_500ctrl);
plot(original, 'LineWidth', 1.2);
legend("Spline approximated signal", "Original signal");
title("a) p=3; No.CP=500");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(3)
hold on;
plot(fn_form_200ctrl);
plot(original, 'LineWidth', 1.2);
legend("Spline approximated signal", "Original signal");
title("b) p=3; No.CP=200");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(4)
hold on;
plot(fn_form_100ctrl);
plot(original, 'LineWidth', 1.2);
legend("Spline approximated signal", "Original signal");
title("c) p=3; No.CP=100");
xlabel("Samples");
ylabel("Magnitude");
hold off;

figure(5)
hold on;
plot(fn_form_50ctrl);
plot(original, 'LineWidth', 1.2);
legend("Spline approximated signal", "Original signal");
title("d) p=3; No.CP=50");
xlabel("Samples");
ylabel("Magnitude");
hold off;

%Performance metrics calculation

signal_vector = [noisy_signal' fn_form_500ctrl' fn_form_200ctrl' fn_form_100ctrl' fn_form_50ctrl']; %From this point onwards this order is used for the signals inside the vectors.

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

t = 0:1:596;
Mocap = amc_to_matrix("11_01.amc");

%Extraction of the elements used on the original paper
ch1 = Mocap(:, 2); %Called channel 1 in the original paper
ch1 = (6*(ch1-17.5)+86)'; %Replicating scaling and offset applied in the original paper
ch3 = Mocap(:, 6)'; %Called channel 3 in the original paper

% ch1 B-Spline filtering

%Spline function fitting using the spap2 function
spline_app_598ctrl_ch1 = spap2(598, 3, t, ch1);
spline_app_200ctrl_ch1 = spap2(200, 3, t, ch1);
spline_app_100ctrl_ch1 = spap2(100, 3, t, ch1);
spline_app_50ctrl_ch1 = spap2(50, 3, t, ch1);
spline_app_25ctrl_ch1 = spap2(25, 3, t, ch1);

fn_form_598ctrl_ch1 = fnval(fn2fm(spline_app_598ctrl_ch1, 'ppform'), t);
fn_form_200ctrl_ch1 = fnval(fn2fm(spline_app_200ctrl_ch1, 'ppform'), t);
fn_form_100ctrl_ch1 = fnval(fn2fm(spline_app_100ctrl_ch1, 'ppform'), t);
fn_form_50ctrl_ch1 = fnval(fn2fm(spline_app_50ctrl_ch1, 'ppform'), t);
fn_form_25ctrl_ch1 = fnval(fn2fm(spline_app_25ctrl_ch1, 'ppform'), t);

figure(1)
subplot(2, 2, 1);
hold on;
plot(fn_form_598ctrl_ch1, 'LineWidth', 1.2);
plot(ch1, "k");
legend("Spline approximated signal", "Original signal");
title("a) p=3; No.CP=598");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 2);
hold on;
plot(fn_form_100ctrl_ch1, 'LineWidth', 1.2);
plot(ch1, "k");
legend("Spline approximated signal", "Original signal");
title("c) p=3; No.CP=100");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 3);
hold on;
plot(fn_form_50ctrl_ch1, 'LineWidth', 1.2);
plot(ch1, "k");
legend("Spline approximated signal", "Original signal");
title("d) p=3; No.CP=50");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 4);
hold on;
plot(fn_form_25ctrl_ch1, 'LineWidth', 1.2);
plot(ch1, "k");
legend("Spline approximated signal", "Original signal");
title("b) p=3; No.CP=25");
xlabel("Samples");
ylabel("Magnitude");
hold off;

%Performance metrics calculation

signal_vector = [ch1' fn_form_598ctrl_ch1' fn_form_200ctrl_ch1' fn_form_100ctrl_ch1' fn_form_50ctrl_ch1' fn_form_25ctrl_ch1']; %From this point onwards this order is used for the signals inside the vectors.

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

% ch3 B-Spline filtering

%Spline function fitting using the spap2 function
spline_app_598ctrl_ch3 = spap2(598, 3, t, ch3);
spline_app_200ctrl_ch3 = spap2(200, 3, t, ch3);
spline_app_100ctrl_ch3 = spap2(100, 3, t, ch3);
spline_app_50ctrl_ch3 = spap2(50, 3, t, ch3);
spline_app_25ctrl_ch3 = spap2(25, 3, t, ch3);

fn_form_598ctrl_ch3 = fnval(fn2fm(spline_app_598ctrl_ch3, 'ppform'), t);
fn_form_200ctrl_ch3 = fnval(fn2fm(spline_app_200ctrl_ch3, 'ppform'), t);
fn_form_100ctrl_ch3 = fnval(fn2fm(spline_app_100ctrl_ch3, 'ppform'), t);
fn_form_50ctrl_ch3 = fnval(fn2fm(spline_app_50ctrl_ch3, 'ppform'), t);
fn_form_25ctrl_ch3 = fnval(fn2fm(spline_app_25ctrl_ch3, 'ppform'), t);

figure(2)
subplot(2, 2, 1);
hold on;
plot(fn_form_598ctrl_ch3, 'LineWidth', 1.2);
plot(ch3, "k");
legend("Spline approximated signal", "Original signal");
title("a) p=3; No.CP=598");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 2);
hold on;
plot(fn_form_100ctrl_ch3, 'LineWidth', 1.2);
plot(ch3, "k");
legend("Spline approximated signal", "Original signal");
title("c) p=3; No.CP=100");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 3);
hold on;
plot(fn_form_50ctrl_ch3, 'LineWidth', 1.2);
plot(ch3, "k");
legend("Spline approximated signal", "Original signal");
title("d) p=3; No.CP=50");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 4);
hold on;
plot(fn_form_25ctrl_ch3, 'LineWidth', 1.2);
plot(ch3, "k");
legend("Spline approximated signal", "Original signal");
title("b) p=3; No.CP=25");
xlabel("Samples");
ylabel("Magnitude");
hold off;

%Performance metrics calculation

signal_vector = [ch3' fn_form_598ctrl_ch3' fn_form_200ctrl_ch3' fn_form_100ctrl_ch3' fn_form_50ctrl_ch3' fn_form_25ctrl_ch3']; %From this point onwards this order is used for the signals inside the vectors.

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