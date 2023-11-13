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
subplot(2, 2, 1);
hold on;
plot(fn_form_500ctrl);
plot(original, 'LineWidth', 1.2);
legend("Spline approximated signal", "Original signal");
title("a) p=3; No.CP=500");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 2);
hold on;
plot(fn_form_200ctrl);
plot(original, 'LineWidth', 1.2);
legend("Spline approximated signal", "Original signal");
title("b) p=3; No.CP=200");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 3);
hold on;
plot(fn_form_100ctrl);
plot(original, 'LineWidth', 1.2);
legend("Spline approximated signal", "Original signal");
title("c) p=3; No.CP=100");
xlabel("Samples");
ylabel("Magnitude");
hold off;

subplot(2, 2, 4);
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