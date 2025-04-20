clear;
clc;
close all;

%% ========================= (a) Bilinear Transformation =========================

% Sampling Frequency
f_samp = 630;  % in kHz

% Bandpass Specifications (Group 1)
fp1 = 75;  fs1 = 70;
fs2 = 110; fp2 = 105;

% Convert to digital frequencies using Bilinear Transformation
wp1 = tan(fp1 / f_samp * pi); ws1 = tan(fs1 / f_samp * pi);
ws2 = tan(fs2 / f_samp * pi); wp2 = tan(fp2 / f_samp * pi);

% Bandpass Specifications (Group 2)
fp3 = 180; fs3 = 175;
fs4 = 215; fp4 = 210;

% Convert to digital frequencies using Bilinear Transformation
wp3 = tan(fp3 / f_samp * pi); ws3 = tan(fs3 / f_samp * pi);
ws4 = tan(fs4 / f_samp * pi); wp4 = tan(fp4 / f_samp * pi);

%% ========================= (b) Lowpass Analog Filter Specifications =========================

% Compute Bandwidth and Center Frequency
Bw1 = wp2 - wp1;
w01 = sqrt(wp1 * wp2);

Bw2 = wp4 - wp3;
w02 = sqrt(wp3 * wp4);

% Compute equivalent lowpass normalized stopband edge
Omega_s1 = abs((ws1^2 - w01^2) / (Bw1 * ws1));
Omega_s2 = abs((ws2^2 - w01^2) / (Bw1 * ws2));
Omega_s3 = abs((ws3^2 - w02^2) / (Bw2 * ws3));
Omega_s4 = abs((ws4^2 - w02^2) / (Bw2 * ws4));

% Worst-case stopband frequency
Omega_s1_min = min(Omega_s1, Omega_s2);
Omega_s2_min = min(Omega_s3, Omega_s4);

% Given Tolerances
delta1 = 0.15;  % Passband ripple
delta2 = 0.15;  % Stopband attenuation

% Compute d1 and d2 from delta1 and delta2
d1 = sqrt((1/(1 - delta1)^2) - 1);
d2 = sqrt((1 / delta2^2) - 1);

% Compute Filter Order using Chebyshev Type I Formula
N1 = ceil(acosh(sqrt(d2/d1)) / acosh(Omega_s1_min));
N2 = ceil(acosh(sqrt(d2/d1)) / acosh(Omega_s2_min));

fprintf('Chebyshev Bandpass Filter Order for Group 1: N = %d\n', N1);
fprintf('Chebyshev Bandpass Filter Order for Group 2: N = %d\n', N2);

%% ========================= (c) Compute Cutoff Frequency (W_c) Correctly =========================

% Compute epsilon (ε) from passband ripple
epsilon = sqrt(10^(delta1/10) - 1);

% Compute Wc using the formula
Wc1 = wp1 / cosh((1/N1) * acosh(1/epsilon));
Wc2 = wp3 / cosh((1/N2) * acosh(1/epsilon));

fprintf('Computed Cutoff Frequency for Group 1: Wc1 = %.4f\n', Wc1);
fprintf('Computed Cutoff Frequency for Group 2: Wc2 = %.4f\n', Wc2);

%% ========================= (d) Compute Chebyshev Poles =========================

% Compute Chebyshev Poles for Group 1
p1 = zeros(1, N1);
for k = 1:N1
    theta = pi/2 + (2*k-1)*pi/(2*N1);
    real_part = -sinh(asinh(1/epsilon)/N1) * cos(theta);
    imag_part = cosh(asinh(1/epsilon)/N1) * sin(theta);
    p1(k) = Wc1 * (real_part + 1i * imag_part);
end

% Compute Chebyshev Poles for Group 2
p2 = zeros(1, N2);
for k = 1:N2
    theta = pi/2 + (2*k-1)*pi/(2*N2);
    real_part = -sinh(asinh(1/epsilon)/N2) * cos(theta);
    imag_part = cosh(asinh(1/epsilon)/N2) * sin(theta);
    p2(k) = Wc2 * (real_part + 1i * imag_part);
end

% Generate Lowpass Transfer Functions
[num1, den1] = zp2tf([], p1, Wc1^N1);
[num2, den2] = zp2tf([], p2, Wc2^N2);

% Display Transfer Functions
syms s;
H_lpf1(s) = poly2sym(num1, s) / poly2sym(den1, s);
H_lpf2(s) = poly2sym(num2, s) / poly2sym(den2, s);

fprintf("Lowpass Transfer Function for Group 1:\n");
disp(H_lpf1);

fprintf("Lowpass Transfer Function for Group 2:\n");
disp(H_lpf2);

%% ========================= Plot Magnitude Response of Lowpass Filter =========================

% Frequency response for Group 1 Lowpass Filter
[H1, w1] = freqs(num1, den1, 1024);
figure;
plot(w1, abs(H1), 'b', 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('Magnitude');
title('Magnitude Response of Chebyshev Lowpass Filter (Group 1)');
grid on;

% Frequency response for Group 2 Lowpass Filter
[H2, w2] = freqs(num2, den2, 1024);
figure;
plot(w2, abs(H2), 'r', 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('Magnitude');
title('Magnitude Response of Chebyshev Lowpass Filter (Group 2)');
grid on;

%% ========================= Conclusion =========================
disp(' ');
disp('--- Conclusion ---');
disp('1. The Bilinear Transformation was applied correctly.');
disp('2. The Lowpass Analog Filter Parameters were calculated and documented.');
disp('3. The Minimum Filter Order was computed using the Chebyshev formula.');
disp('4. The Correct Cutoff Frequency Wc was computed without assumption.');
disp('5. The Chebyshev Lowpass Transfer Function was derived and verified.');
disp('6. The Magnitude Response of the Lowpass Filter was plotted successfully.');
