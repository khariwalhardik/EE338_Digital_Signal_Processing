% Sampling Frequency
f_samp = 630;

% Bandpass Specifications
fp1 = 180; fs1 = 175;
fs2 = 215; fp2 = 210;

% Convert to digital frequencies
wp1 = tan(fp1/f_samp*pi); ws1 = tan(fs1/f_samp*pi);
ws2 = tan(fs2/f_samp*pi); wp2 = tan(fp2/f_samp*pi);

% Butterworth Bandpass Design Parameters
Wc = 1.08; N = 10;

% Analog LPF Poles for Butterworth Filter
p = zeros(1, N);  
for k = 1:N
    theta = pi/2 + (2*k-1)*pi/(2*N);
    p(k) = Wc * (cos(theta) + 1i*sin(theta));
end
[num, den] = zp2tf([], p, Wc^N);

syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);

% Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;
analog_bpf(s) = analog_lpf((s*s + W0*W0)/(B*s)); 
discrete_bpf(z) = analog_bpf((z-1)/(z+1));

% Coefficients of Discrete BPF
[nz_bpf, dz_bpf] = numden(discrete_bpf(z));                   
nz_bpf = sym2poly(expand(nz_bpf));
dz_bpf = sym2poly(expand(dz_bpf));                            
k_bpf = dz_bpf(1);                                            
dz_bpf = dz_bpf/k_bpf;
nz_bpf = nz_bpf/k_bpf;

% Frequency Response of Bandpass Filter
[H_bpf, w_bpf] = freqz(nz_bpf, dz_bpf, 1024, f_samp);
figure;
plot(w_bpf, abs(H_bpf));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Response of Bandpass Filter (Group-II)');
grid on;

% Frequencies of Interest
freqs = [175, 180, 210, 215];

% Find Indices Closest to Desired Frequencies
[~, idx] = arrayfun(@(f) min(abs(w_bpf - f)), freqs);


% Display Frequency Values
disp('Magnitude Values at Specific Frequencies:');
for i = 1:length(freqs)
    fprintf('Frequency %d Hz: Magnitude = %.4f\n', freqs(i), abs(H_bpf(idx(i))));
end
