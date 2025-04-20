clear;
clc;
close all;

%% Group-II Bandpass Filter Design
f_samp_g2 = 630;  % Sampling Frequency

% Bandpass Specifications for Group-II
fp1_g2 = 180; fs1_g2 = 175;
fs2_g2 = 215; fp2_g2 = 210;

% Convert to digital frequencies for Group-II
wp1_g2 = tan(fp1_g2/f_samp_g2*pi); ws1_g2 = tan(fs1_g2/f_samp_g2*pi);
ws2_g2 = tan(fs2_g2/f_samp_g2*pi); wp2_g2 = tan(fp2_g2/f_samp_g2*pi);

% Butterworth Bandpass Design Parameters for Group-II
Wc_g2 = 1.08; N_g2 = 11;

% Analog LPF Poles for Butterworth Filter (Group-II)
p_g2 = zeros(1, N_g2);  
for k = 1:N_g2
    theta_g2 = pi/2 + (2*k-1)*pi/(2*N_g2);
    p_g2(k) = Wc_g2 * (cos(theta_g2) + 1i*sin(theta_g2));
end
[num_g2, den_g2] = zp2tf([], p_g2, Wc_g2^N_g2);

syms s z;
analog_lpf_g2(s) = poly2sym(num_g2, s) / poly2sym(den_g2, s);

% Bandpass Transformation for Group-II
W0_g2 = sqrt(wp1_g2 * wp2_g2);
B_g2 = wp2_g2 - wp1_g2;
analog_bpf_g2(s) = analog_lpf_g2((s^2 + W0_g2^2) / (B_g2 * s)); 
discrete_bpf_g2(z) = analog_bpf_g2((z-1) / (z+1));

% Coefficients of Discrete BPF for Group-II
[nz_bpf_g2, dz_bpf_g2] = numden(discrete_bpf_g2(z));                   
nz_bpf_g2 = sym2poly(expand(nz_bpf_g2));
dz_bpf_g2 = sym2poly(expand(dz_bpf_g2));                            
k_bpf_g2 = dz_bpf_g2(1);                                            
dz_bpf_g2 = dz_bpf_g2 / k_bpf_g2;
nz_bpf_g2 = nz_bpf_g2 / k_bpf_g2;

% Frequency Response of Group-II Bandpass Filter
[~, ~] = freqz(nz_bpf_g2, dz_bpf_g2, 1024, f_samp_g2);


%% Group-I Bandpass Filter Design
f_samp_g1 = 630;  % Same Sampling Frequency for Combining

% Bandpass Specifications for Group-I
fp1_g1 = 75; fs1_g1 = 70;
fs2_g1 = 110; fp2_g1 = 105;

% Convert to digital frequencies for Group-I
wp1_g1 = tan(fp1_g1/f_samp_g1*pi); ws1_g1 = tan(fs1_g1/f_samp_g1*pi);
ws2_g1 = tan(fs2_g1/f_samp_g1*pi); wp2_g1 = tan(fp2_g1/f_samp_g1*pi);

% Butterworth Bandpass Design Parameters for Group-I
Wc_g1 = 1.08; N_g1 = 11;

% Analog LPF Poles for Butterworth Filter (Group-I)
p_g1 = zeros(1, N_g1);  
for k = 1:N_g1
    theta_g1 = pi/2 + (2*k-1)*pi/(2*N_g1);
    p_g1(k) = Wc_g1 * (cos(theta_g1) + 1i*sin(theta_g1));
end
[num_g1, den_g1] = zp2tf([], p_g1, Wc_g1^N_g1);

analog_lpf_g1(s) = poly2sym(num_g1, s) / poly2sym(den_g1, s);

% Bandpass Transformation for Group-I
W0_g1 = sqrt(wp1_g1 * wp2_g1);
B_g1 = wp2_g1 - wp1_g1;
analog_bpf_g1(s) = analog_lpf_g1((s^2 + W0_g1^2) / (B_g1 * s)); 
discrete_bpf_g1(z) = analog_bpf_g1((z-1) / (z+1));

% Coefficients of Discrete BPF for Group-I
[nz_bpf_g1, dz_bpf_g1] = numden(discrete_bpf_g1(z));                   
nz_bpf_g1 = sym2poly(expand(nz_bpf_g1));
dz_bpf_g1 = sym2poly(expand(dz_bpf_g1));                            
k_bpf_g1 = dz_bpf_g1(1);                                            
dz_bpf_g1 = dz_bpf_g1 / k_bpf_g1;
nz_bpf_g1 = nz_bpf_g1 / k_bpf_g1;

% Calculate Frequency Response for Both Bandpass Filters
[H_bpf_g1, w_bpf_g1] = freqz(nz_bpf_g1, dz_bpf_g1, 4096, f_samp_g1);
[H_bpf_g2, w_bpf_g2] = freqz(nz_bpf_g2, dz_bpf_g2, 4096, f_samp_g2);

% Ensure Frequency Vectors are the Same
w_combined = w_bpf_g1; % They should be the same if f_samp is identical

% Normalize the Frequency Responses
H_bpf_g1 = abs(H_bpf_g1) / max(abs(H_bpf_g1));
H_bpf_g2 = abs(H_bpf_g2) / max(abs(H_bpf_g2));

% Combine the Two Bandpass Filters (Multipass Band Filter)
H_multipass = max(H_bpf_g1, H_bpf_g2); % Use max to combine

% Plot for Group-I BPF
figure;
plot(w_combined, H_bpf_g1, 'b', 'LineWidth', 1.5);
title('Group-I Bandpass Filter Response');
xlabel('Frequency (KHz)');
ylabel('Magnitude');
grid on;

% Plot for Group-II BPF
figure;
plot(w_combined, H_bpf_g2, 'g', 'LineWidth', 1.5);
title('Group-II Bandpass Filter Response');
xlabel('Frequency (KHz)');
ylabel('Magnitude');
grid on;

% Plot for Combined Multipass BPF
figure;
plot(w_combined, H_multipass, 'r', 'LineWidth', 1.5);
title('Multipass Band Filter Response');
xlabel('Frequency (KHz)');
ylabel('Magnitude');
grid on;

% Group-I Passband and Stopband Frequencies
frequencies_g1 = [fs1_g1, fp1_g1, fp2_g1, fs2_g1];
magnitudes_g1 = interp1(w_combined, H_bpf_g1, frequencies_g1);

% Group-II Passband and Stopband Frequencies
frequencies_g2 = [fs1_g2, fp1_g2, fp2_g2, fs2_g2];
magnitudes_g2 = interp1(w_combined, H_bpf_g2, frequencies_g2);

% Display Magnitudes for Group-I
disp('Group-I Passband and Stopband Magnitudes:');
for i = 1:length(frequencies_g1)
    fprintf('Frequency = %d KHz, Magnitude = %.4f\n', frequencies_g1(i), magnitudes_g1(i));
end

% Display Magnitudes for Group-II
disp('Group-II Passband and Stopband Magnitudes:');
for i = 1:length(frequencies_g2)
    fprintf('Frequency = %d KHz, Magnitude = %.4f\n', frequencies_g2(i), magnitudes_g2(i));
end
