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
Wc2 = 0.85;  % For Group 2

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

%====================== Print Lowpass Filter Transfer Functions ============================
disp('Analog Lowpass Filter Transfer Function (Group 1):');
pretty(vpa(analog_lpf1(s), 4));

disp('Analog Lowpass Filter Transfer Function (Group 2):');
pretty(vpa(analog_lpf2(s), 4));

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

% Transform Lowpass to Bandpass (Group 2)
analog_bpf2(s) = analog_lpf2((s^2 + W0_2^2) / (B2 * s));

%====================== Print Bandpass Filter Transfer Functions ============================
disp('Analog Bandpass Filter Transfer Function (Group 1):');
pretty(vpa(analog_bpf1(s), 4));

disp('Analog Bandpass Filter Transfer Function (Group 2):');
pretty(vpa(analog_bpf2(s), 4));

%====================== Bilinear Transformation to Discrete Domain ============================
syms z;

% Bilinear Transform: s = (2/T) * (z-1)/(z+1)
T = 1 / f_samp;  % Sampling period
bilinear_s = (2/T) * (z - 1) / (z + 1);

% Apply Bilinear Transformation
discrete_bpf1(z) = simplify(subs(analog_bpf1(s), s, bilinear_s));
discrete_bpf2(z) = simplify(subs(analog_bpf2(s), s, bilinear_s));

%====================== Print Discrete-Time Transfer Functions ============================
disp('Discrete-Time Bandpass Filter Transfer Function (Group 1):');
pretty(vpa(discrete_bpf1(z), 4));

disp('Discrete-Time Bandpass Filter Transfer Function (Group 2):');
pretty(vpa(discrete_bpf2(z), 4));

%====================== Discrete Frequency Response ============================
w_dt = linspace(0, pi, 1000);  % Discrete frequency range

% Check for division by zero in the denominator
[num1_z, den1_z] = numden(discrete_bpf1(z));
[num2_z, den2_z] = numden(discrete_bpf2(z));

% Compute denominator values to check for zero
den1_eval = double(subs(den1_z, z, exp(1i * w_dt)));
den2_eval = double(subs(den2_z, z, exp(1i * w_dt)));

% Identify problematic frequencies
zero_idx1 = find(abs(den1_eval) < 1e-6);
zero_idx2 = find(abs(den2_eval) < 1e-6);

if ~isempty(zero_idx1)
    disp('Warning: Division by zero detected in Group 1 at frequencies:');
    disp(w_dt(zero_idx1));
end
if ~isempty(zero_idx2)
    disp('Warning: Division by zero detected in Group 2 at frequencies:');
    disp(w_dt(zero_idx2));
end

% Avoid division by zero by adding a small epsilon
epsilon = 1e-3;
disp('Discrete BPF1:');
disp(discrete_bpf1(z));
disp('w_dt values:');
disp(w_dt);

H_dt1 = abs(double(subs(discrete_bpf1(z), z, exp(1i * w_dt)) + epsilon));
H_dt2 = abs(double(subs(discrete_bpf2(z), z, exp(1i * w_dt)) + epsilon));

% Plot Discrete-Time Frequency Response
figure;
plot(w_dt, H_dt1, 'r', 'LineWidth', 1.5);
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
title('Magnitude Response of Discrete Bandpass Filter (Group 1)');
grid on;

figure;
plot(w_dt, H_dt2, 'm', 'LineWidth', 1.5);
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
title('Magnitude Response of Discrete Bandpass Filter (Group 2)');
grid on;
