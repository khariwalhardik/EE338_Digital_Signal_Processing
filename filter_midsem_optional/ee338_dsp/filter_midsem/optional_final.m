clear;
clc;
close all;

% Sampling Frequency
f_samp = 630;

%====================== Group 1 Bandpass Filter ============================
% Bandpass Specifications
fp1 = 75; fs1 = 70;
fs2 = 110; fp2 = 105;

% Convert to digital frequencies using Bilinear Transformation
wp1 = tan(fp1 / f_samp * pi); ws1 = tan(fs1 / f_samp * pi);
ws2 = tan(fs2 / f_samp * pi); wp2 = tan(fp2 / f_samp * pi);

% Chebyshev Type I Filter Design Parameters
Wc1 = 0.94; N1 = 5; Rp = 0.15;  % Passband ripple in dB
epsilon = sqrt(10^(Rp/10) - 1); % Compute epsilon from passband ripple

% Compute Chebyshev Poles
p1 = zeros(1, N1);
for k = 1:N1
    theta = pi/2 + (2*k-1)*pi/(2*N1);
    real_part = -sinh(asinh(1/epsilon)/N1) * cos(theta);
    imag_part = cosh(asinh(1/epsilon)/N1) * sin(theta);
    p1(k) = Wc1 * (real_part + 1i * imag_part);
end

[num1, den1] = zp2tf([], p1, Wc1^N1);
syms s z;
analog_lpf1(s) = poly2sym(num1, s) / poly2sym(den1, s);

% Bandpass Transformation
W0_1 = sqrt(wp1 * wp2);
B1 = wp2 - wp1;
analog_bpf1(s) = analog_lpf1((s^2 + W0_1^2) / (B1 * s));
discrete_bpf1(z) = analog_bpf1((z - 1) / (z + 1));

% Coefficients of Discrete BPF
[nz_bpf1, dz_bpf1] = numden(discrete_bpf1(z));                   
nz_bpf1 = sym2poly(expand(nz_bpf1));
dz_bpf1 = sym2poly(expand(dz_bpf1));                            
k_bpf1 = dz_bpf1(1);                                            
dz_bpf1 = dz_bpf1 / k_bpf1;
nz_bpf1 = nz_bpf1 / k_bpf1;

% Frequency Response of Bandpass Filter
[H_bpf1, w_bpf] = freqz(nz_bpf1, dz_bpf1, 1024, f_samp);
H_bpf1 = H_bpf1 / max(abs(H_bpf1)); % Normalize magnitude to keep it ≤ 1

%====================== Group 2 Bandpass Filter ============================
% Bandpass Specifications
fp1 = 180; fs1 = 175;
fs2 = 215; fp2 = 210;

% Convert to digital frequencies using Bilinear Transformation
wp1 = tan(fp1 / f_samp * pi); ws1 = tan(fs1 / f_samp * pi);
ws2 = tan(fs2 / f_samp * pi); wp2 = tan(fp2 / f_samp * pi);

% Chebyshev Type I Filter Design Parameters
Wc2 = 0.95; N2 = 5; Rp = 0.15;  % Passband ripple in dB

% Compute Chebyshev Poles
p2 = zeros(1, N2);
for k = 1:N2
    theta = pi/2 + (2*k-1)*pi/(2*N2);
    real_part = -sinh(asinh(1/epsilon)/N2) * cos(theta);
    imag_part = cosh(asinh(1/epsilon)/N2) * sin(theta);
    p2(k) = Wc2 * (real_part + 1i * imag_part);
end

[num2, den2] = zp2tf([], p2, Wc2^N2);
analog_lpf2(s) = poly2sym(num2, s) / poly2sym(den2, s);

% Bandpass Transformation
W0_2 = sqrt(wp1 * wp2);
B2 = wp2 - wp1;
analog_bpf2(s) = analog_lpf2((s^2 + W0_2^2) / (B2 * s));
discrete_bpf2(z) = analog_bpf2((z - 1) / (z + 1));

% Coefficients of Discrete BPF
[nz_bpf2, dz_bpf2] = numden(discrete_bpf2(z));                   
nz_bpf2 = sym2poly(expand(nz_bpf2));
dz_bpf2 = sym2poly(expand(dz_bpf2));                            
k_bpf2 = dz_bpf2(1);                                            
dz_bpf2 = dz_bpf2 / k_bpf2;
nz_bpf2 = nz_bpf2 / k_bpf2;

% Frequency Response of Bandpass Filter
[H_bpf2, ~] = freqz(nz_bpf2, dz_bpf2, 1024, f_samp);
H_bpf2 = H_bpf2 / max(abs(H_bpf2)); % Normalize magnitude to keep it ≤ 1

%====================== Multi-Bandpass Filter ============================
H_multi = abs(H_bpf1) + abs(H_bpf2); % Combining both bandpass filters

%====================== Pole-Zero Plot ============================
figure;
subplot(1,3,1);
plot(real(p1), imag(p1), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
circle = exp(1i * linspace(0, 2*pi, 100));
plot(real(circle), imag(circle), 'k--');
hold off;
title('Poles of Group-I Bandpass Filter');
xlabel('Real'); ylabel('Imaginary');
axis equal; grid on;

subplot(1,3,2);
plot(real(p2), imag(p2), 'gx', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
plot(real(circle), imag(circle), 'k--');
hold off;
title('Poles of Group-II Bandpass Filter');
xlabel('Real'); ylabel('Imaginary');
axis equal; grid on;

subplot(1,3,3);
plot(real([p1, p2]), imag([p1, p2]), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
plot(real(circle), imag(circle), 'k--');
hold off;
title('Poles of Combined Multi-Bandpass Filter');
xlabel('Real'); ylabel('Imaginary');
axis equal; grid on;

%====================== Display Magnitude Values ============================
freqs = [70, 75, 105, 110, 175, 180, 210, 215];

% Find Indices Closest to Desired Frequencies
[~, idx] = arrayfun(@(f) min(abs(w_bpf - f)), freqs);

% Display Frequency Values
disp('Magnitude Values at Specific Frequencies:');
for i = 1:length(freqs)
    fprintf('Frequency %d Hz: Magnitude = %.4f\n', freqs(i), H_multi(idx(i)));
end

%====================== Plotting Responses ============================
figure;
plot(w_bpf, abs(H_bpf1), 'b', 'LineWidth', 1.5);
xlabel('Frequency (KHz)'); ylabel('Magnitude');
title('Magnitude Response of Group 1 Bandpass Filter');
grid on;

figure;
plot(w_bpf, abs(H_bpf2), 'g', 'LineWidth', 1.5);
xlabel('Frequency (KHz)'); ylabel('Magnitude');
title('Magnitude Response of Group 2 Bandpass Filter');
grid on;

figure;
plot(w_bpf, H_multi, 'r', 'LineWidth', 2);
xlabel('Frequency (KHz)'); ylabel('Magnitude');
title('Magnitude Response of Multi-Bandpass Filter');
legend('Combined Multi-Bandpass');
grid on;
