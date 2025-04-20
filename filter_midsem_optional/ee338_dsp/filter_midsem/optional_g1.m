clear;
clc;
close all;

% Sampling Frequency
f_samp = 630;

% Bandpass Specifications
fp1 = 75; fs1 = 70;
fs2 = 110; fp2 = 105;

% Convert to digital frequencies using Bilinear Transformation
wp1 = tan(fp1 / f_samp * pi); ws1 = tan(fs1 / f_samp * pi);
ws2 = tan(fs2 / f_samp * pi); wp2 = tan(fp2 / f_samp * pi);

% Chebyshev Type I Filter Design Parameters
Wc = 1; N = 5; Rp = 0.15;  % Passband ripple in dB

epsilon = sqrt(10^(Rp/10) - 1); % Compute epsilon from passband ripple

% Compute Chebyshev Poles
p = zeros(1, N);
for k = 1:N
    theta = pi/2 + (2*k-1)*pi/(2*N);
    real_part = -sinh(asinh(1/epsilon)/N) * cos(theta);
    imag_part = cosh(asinh(1/epsilon)/N) * sin(theta);
    p(k) = Wc * (real_part + 1i * imag_part);
end

[num, den] = zp2tf([], p, Wc^N);

syms s z;
analog_lpf(s) = poly2sym(num, s) / poly2sym(den, s);

% Bandpass Transformation
W0 = sqrt(wp1 * wp2);
B = wp2 - wp1;
analog_bpf(s) = analog_lpf((s^2 + W0^2) / (B * s));
discrete_bpf(z) = analog_bpf((z - 1) / (z + 1));

% Coefficients of Discrete BPF
[nz_bpf, dz_bpf] = numden(discrete_bpf(z));                   
nz_bpf = sym2poly(expand(nz_bpf));
dz_bpf = sym2poly(expand(dz_bpf));                            
k_bpf = dz_bpf(1);                                            
dz_bpf = dz_bpf / k_bpf;
nz_bpf = nz_bpf / k_bpf;

% Frequency Response of Bandpass Filter
[H_bpf, w_bpf] = freqz(nz_bpf, dz_bpf, 1024, f_samp);
H_bpf = H_bpf / max(abs(H_bpf)); % Normalize magnitude to keep it ≤ 1

figure;
plot(w_bpf, abs(H_bpf));
xlabel('Frequency (KHz)');
ylabel('Magnitude');
title('Magnitude Response of Chebyshev Type I Bandpass Filter (Group-I)');
grid on;

% Frequencies of Interest
freqs = [70, 75, 105, 110];

% Find Indices Closest to Desired Frequencies
[~, idx] = arrayfun(@(f) min(abs(w_bpf - f)), freqs);

% Plot Markers at Specified Frequencies
hold on;

% Display Frequency Values
disp('Magnitude Values at Specific Frequencies:');
for i = 1:length(freqs)
    fprintf('Frequency %d Hz: Magnitude = %.4f\n', freqs(i), abs(H_bpf(idx(i))));
end
