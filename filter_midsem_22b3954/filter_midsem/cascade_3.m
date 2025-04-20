% Specifications for Bandpass Filter
fp1_bp = 75e3;
fp2_bp = 210e3;
fs1_bp = 70e3;
fs2_bp = 215e3;

% Specifications for Bandstop Filter
fs1_bs = 105e3;
fs2_bs = 180e3;
fp1_bs = 100e3;
fp2_bs = 185e3;

% Common Parameters
f_samp = 630e3;
Wc = 1.08; 
N = 8;

% Butterworth Poles Calculation
p = zeros(1, N);
for k = 1:N
    theta = pi/2 + (2*k-1)*pi/(2*N);
    p(k) = Wc * (cos(theta) + 1i*sin(theta));
end

% Analog LPF Transfer Function
[num, den] = zp2tf([], p, Wc^N);

% Bandpass Transformation
wp1_bp = tan(fp1_bp/f_samp*pi);
wp2_bp = tan(fp2_bp/f_samp*pi);
W0_bp = sqrt(wp1_bp * wp2_bp);
B_bp = wp2_bp - wp1_bp;
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);
analog_bpf(s) = analog_lpf((s*s + W0_bp^2)/(B_bp*s));
discrete_bpf(z) = analog_bpf((z-1)/(z+1));

% Bandstop Transformation
wp1_bs = tan(fp1_bs/f_samp*pi);
wp2_bs = tan(fp2_bs/f_samp*pi);
W0_bs = sqrt(wp1_bs * wp2_bs);
B_bs = wp2_bs - wp1_bs;
analog_bsf(s) = analog_lpf((B_bs*s)/(s*s + W0_bs^2));
discrete_bsf(z) = analog_bsf((z-1)/(z+1));

% Coefficients of Discrete BPF
[nz_bpf, dz_bpf] = numden(discrete_bpf(z));
nz_bpf = sym2poly(expand(nz_bpf));
dz_bpf = sym2poly(expand(dz_bpf));
k_bpf = dz_bpf(1);
dz_bpf = dz_bpf/k_bpf;
nz_bpf = nz_bpf/k_bpf;

% Coefficients of Discrete BSF
[nz_bsf, dz_bsf] = numden(discrete_bsf(z));
nz_bsf = sym2poly(expand(nz_bsf));
dz_bsf = sym2poly(expand(dz_bsf));
k_bsf = dz_bsf(1);
dz_bsf = dz_bsf/k_bsf;
nz_bsf = nz_bsf/k_bsf;

% Cascading Bandpass and Bandstop Filters
nz_combined = conv(nz_bpf, nz_bsf);
dz_combined = conv(dz_bpf, dz_bsf);

% Frequency Response of Combined Filter
fvtool(nz_combined, dz_combined);

% Magnitude Plot
[H_combined, f] = freqz(nz_combined, dz_combined, 1024*1024, f_samp);
plot(f, abs(H_combined));
title('Magnitude Response of Combined Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
