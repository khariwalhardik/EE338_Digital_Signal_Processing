clear;
clc;
close all;

% Sampling Frequency
f_samp = 630;

%====================== Chebyshev Lowpass Filter Design ============================
% Given Passband Ripple and Stopband Attenuation
Rp = 0.15;  % Passband ripple in dB

% Compute Epsilon for Chebyshev Filter
epsilon = sqrt(10^(Rp/10) - 1);

% Compute Cutoff Frequencies for Both Groups
Wc1 = 1.08;  % For Group 1
Wc2 = 2.50;  % For Group 2

% Filter Order for Both Groups
N1 = 5;  % For Group 1
N2 = 5;  % For Group 2

% Compute Poles of Chebyshev Lowpass Filter (Group 1)
p1 = zeros(1, N1);
for k = 1:N1
    theta = pi/2 + (2*k-1)*pi/(2*N1);
    real_part = -sinh(asinh(1/epsilon)/N1) * cos(theta);
    imag_part = cosh(asinh(1/epsilon)/N1) * sin(theta);
    p1(k) = Wc1 * (real_part + 1i * imag_part);
end

% Compute Poles of Chebyshev Lowpass Filter (Group 2)
p2 = zeros(1, N2);
for k = 1:N2
    theta = pi/2 + (2*k-1)*pi/(2*N2);
    real_part = -sinh(asinh(1/epsilon)/N2) * cos(theta);
    imag_part = cosh(asinh(1/epsilon)/N2) * sin(theta);
    p2(k) = Wc2 * (real_part + 1i * imag_part);
end

% Compute Transfer Function Coefficients
[num1, den1] = zp2tf([], p1, Wc1^N1);
[num2, den2] = zp2tf([], p2, Wc2^N2);
syms s;
analog_lpf1(s) = poly2sym(num1, s) / poly2sym(den1, s);
analog_lpf2(s) = poly2sym(num2, s) / poly2sym(den2, s);

%====================== Plot Analog Lowpass Filter Responses (Separate Figures) ============================
w = linspace(0.01, 3, 1000);
H_lpf1 = abs(double(subs(analog_lpf1(s), s, 1i*w)));
H_lpf2 = abs(double(subs(analog_lpf2(s), s, 1i*w)));

% Plot Lowpass Filter Response (Group 1)
figure;
plot(w, H_lpf1/3, 'b', 'LineWidth', 1.5);
xlabel('Frequency (rad/sec)');
ylabel('Magnitude');
title('Magnitude Response of Chebyshev Lowpass Filter (Group 1)');
grid on;

% Plot Lowpass Filter Response (Group 2)
figure;
plot(w, H_lpf2/3, 'g', 'LineWidth', 1.5);
xlabel('Frequency (rad/sec)');
ylabel('Magnitude');
title('Magnitude Response of Chebyshev Lowpass Filter (Group 2)');
grid on;

%====================== Bandpass Transformation for Both Groups ============================
% Group 1 Bandpass Specifications
fp1_g1 = 75; fs1_g1 = 70;
fs2_g1 = 110; fp2_g1 = 105;

% Convert to Normalized Digital Frequencies
wp1_g1 = tan(fp1_g1 / f_samp * pi); ws1_g1 = tan(fs1_g1 / f_samp * pi);
ws2_g1 = tan(fs2_g1 / f_samp * pi); wp2_g1 = tan(fp2_g1 / f_samp * pi);

% Compute Center Frequency and Bandwidth
W0_1 = sqrt(wp1_g1 * wp2_g1);
B1 = wp2_g1 - wp1_g1;

% Ensure Bandwidth is not zero
if B1 == 0
    error('Bandwidth B1 cannot be zero. Check wp1 and wp2 values.');
end

% Transform Lowpass to Bandpass (Group 1)
analog_bpf1(s) = analog_lpf1((s^2 + W0_1^2) / (B1 * s));

% Group 2 Bandpass Specifications
fp1_g2 = 180; fs1_g2 = 175;
fs2_g2 = 215; fp2_g2 = 210;

% Convert to Normalized Digital Frequencies
wp1_g2 = tan(fp1_g2 / f_samp * pi); ws1_g2 = tan(fs1_g2 / f_samp * pi);
ws2_g2 = tan(fs2_g2 / f_samp * pi); wp2_g2 = tan(fp2_g2 / f_samp * pi);

% Compute Center Frequency and Bandwidth
W0_2 = sqrt(wp1_g2 * wp2_g2);
B2 = wp2_g2 - wp1_g2;

% Ensure Bandwidth is not zero
if B2 == 0
    error('Bandwidth B2 cannot be zero. Check wp1 and wp2 values.');
end

% Transform Lowpass to Bandpass (Group 2)
analog_bpf2(s) = analog_lpf2((s^2 + W0_2^2) / (B2 * s));

%====================== Plot Analog Bandpass Filter Responses (Separate Figures) ============================
w_bp = linspace(0.01, 3, 1000);
H_bpf1 = abs(double(subs(analog_bpf1(s), s, 1i*w_bp)));
H_bpf2 = abs(double(subs(analog_bpf2(s), s, 1i*w_bp)));

% Plot Bandpass Filter Response (Group 1)
figure;
plot(w_bp, H_bpf1/3, 'r', 'LineWidth', 1.5);
xlabel('Frequency (rad/sec)');
ylabel('Magnitude');
title('Magnitude Response of Chebyshev Bandpass Filter (Group 1)');
grid on;

% Plot Bandpass Filter Response (Group 2)
figure;
plot(w_bp, H_bpf2/3, 'm', 'LineWidth', 1.5);
xlabel('Frequency (rad/sec)');
ylabel('Magnitude');
title('Magnitude Response of Chebyshev Bandpass Filter (Group 2)');
grid on;
