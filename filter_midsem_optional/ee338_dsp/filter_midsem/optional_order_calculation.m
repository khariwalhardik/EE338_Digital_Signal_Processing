clear;
clc;
close all;

%==============================================================================================
% Sampling Frequency (in kHz)
f_samp = 630;  

% Given Passband and Stopband Edge Frequencies (in kHz) for Group 1
fs1 = 70;   % Lower stopband edge
fp1 = 75;   % Lower passband edge
fp2 = 105;  % Upper passband edge
fs2 = 110;  % Upper stopband edge

% Convert to normalized digital frequencies using Bilinear Transformation
ws1 = tan((fs1 / f_samp) * pi);
wp1 = tan((fp1 / f_samp) * pi);
wp2 = tan((fp2 / f_samp) * pi);
ws2 = tan((fs2 / f_samp) * pi);

% Compute Bandwidth and Center Frequency
Bw = wp2 - wp1;  % Bandwidth
w0 = sqrt(wp1 * wp2);  % Center frequency

% Compute equivalent lowpass normalized stopband edge
Omega_s1 = (ws1^2 - w0^2) / (Bw * ws1);
Omega_s2 = (ws2^2 - w0^2) / (Bw * ws2);

% Choose the worst-case normalized stopband frequency
Omega_s = min(abs(Omega_s1), abs(Omega_s2));

% Given Tolerances
delta1 = 0.15;  % Passband ripple
delta2 = 0.15;  % Stopband attenuation

% Compute d1 and d2 from delta1 and delta2
d1 = sqrt((1/(1 - delta1)^2) - 1);
d2 = sqrt((1 / delta2^2) - 1);

% Compute Chebyshev Type I Filter Order using cosh inverse formula
N_chebyshev = ceil(acosh(sqrt(d2/d1)) / acosh(Omega_s));

% Compute Wc (Cutoff Frequency for Lowpass Prototype)
Wc = w0 / cosh((1/N_chebyshev) * acosh(Omega_s));

% Display Filter Order and Wc
fprintf('Chebyshev Bandpass Filter Order Group 1: N = %d\n', N_chebyshev);
fprintf('Cutoff Frequency Wc Group 1: %.4f\n', Wc);

%==============================================================================================
% Given Passband and Stopband Edge Frequencies (in kHz) for Group 2
fs1 = 175;   % Lower stopband edge
fp1 = 180;   % Lower passband edge
fp2 = 210;  % Upper passband edge
fs2 = 215;  % Upper stopband edge

% Convert to normalized digital frequencies using Bilinear Transformation
ws1 = tan((fs1 / f_samp) * pi);
wp1 = tan((fp1 / f_samp) * pi);
wp2 = tan((fp2 / f_samp) * pi);
ws2 = tan((fs2 / f_samp) * pi);

% Compute Bandwidth and Center Frequency
Bw = wp2 - wp1;  % Bandwidth
w0 = sqrt(wp1 * wp2);  % Center frequency

% Compute equivalent lowpass normalized stopband edge
Omega_s1 = (ws1^2 - w0^2) / (Bw * ws1);
Omega_s2 = (ws2^2 - w0^2) / (Bw * ws2);

% Choose the worst-case normalized stopband frequency
Omega_s = min(abs(Omega_s1), abs(Omega_s2));

% Compute Chebyshev Type I Filter Order using cosh inverse formula
N_chebyshev = ceil(acosh(sqrt(d2/d1)) / acosh(Omega_s));

% Compute Wc (Cutoff Frequency for Lowpass Prototype)
Wc = w0 / cosh((1/N_chebyshev) * acosh(Omega_s));

% Display Filter Order and Wc
fprintf('Chebyshev Bandpass Filter Order Group 2: N = %d\n', N_chebyshev);
fprintf('Cutoff Frequency Wc Group 2: %.4f\n', Wc);
