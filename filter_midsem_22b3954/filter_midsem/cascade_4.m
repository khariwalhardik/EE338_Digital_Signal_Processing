% Sampling Frequency
f_samp = 630;

%% Butterworth Analog LPF parameters for Bandpass
Wc_bp = 1.08;   % Cut-off frequency for Bandpass
N_bp = 17;      % Order for Bandpass 

% Calculate poles for Bandpass
p_bp = zeros(1, N_bp); 
for k = 1:N_bp
    theta = pi/2 + (2*k-1)*pi/(2*N_bp);
    p_bp(k) = Wc_bp * (cos(theta) + 1i*sin(theta));
end

%% Bandpass Filter Design
% Bandpass Specifications
fp1 = 75;
fs1 = 70;
fs2 = 215;
fp2 = 210;

% Transformed Band Edge specs using Bilinear Transformation
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

% Parameters for Bandpass Transformation
W0_bp = sqrt(wp1*wp2);
B_bp = wp2-wp1;

% Analog LPF Transfer Function for Bandpass
[num_bp, den_bp] = zp2tf([], p_bp, Wc_bp^N_bp);

% Frequency Response of Bandpass Filter
syms s z;
analog_lpf_bp(s) = poly2sym(num_bp,s)/poly2sym(den_bp,s);      
analog_bpf(s) = analog_lpf_bp((s*s + W0_bp*W0_bp)/(B_bp*s));     
discrete_bpf(z) = analog_bpf((z-1)/(z+1));            

% Coeffs of discrete BPF
[nz_bpf, dz_bpf] = numden(discrete_bpf(z));                   
nz_bpf = sym2poly(expand(nz_bpf));
dz_bpf = sym2poly(expand(dz_bpf));                            
k = dz_bpf(1);                                            
dz_bpf = dz_bpf/k;
nz_bpf = nz_bpf/k;

%% Butterworth Analog LPF parameters for Bandstop
Wc_bs = 1.08;   % Cut-off frequency for Bandstop
N_bs = 16;      % Order for Bandstop 

% Calculate poles for Bandstop
p_bs = zeros(1, N_bs); 
for k = 1:N_bs
    theta = pi/2 + (2*k-1)*pi/(2*N_bs);
    p_bs(k) = Wc_bs * (cos(theta) + 1i*sin(theta));
end

%% Bandstop Filter Design
% Bandstop Specifications
fp1 = 105;
fs1 = 110;
fs2 = 175;
fp2 = 180;

% Transformed Band Edge specs using Bilinear Transformation
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

% Parameters for Bandstop Transformation
W0_bs = sqrt(wp1*wp2);
B_bs = wp2-wp1;

% Analog LPF Transfer Function for Bandstop
[num_bs, den_bs] = zp2tf([], p_bs, Wc_bs^N_bs);

% Frequency Response of Bandstop Filter
analog_lpf_bs(s) = poly2sym(num_bs,s)/poly2sym(den_bs,s);      
analog_bsf(s) = analog_lpf_bs((B_bs*s)/(s*s + W0_bs*W0_bs));     
discrete_bsf(z) = analog_bsf((z-1)/(z+1));            

% Coeffs of discrete BSF
[nz_bsf, dz_bsf] = numden(discrete_bsf(z));                   
nz_bsf = sym2poly(expand(nz_bsf));
dz_bsf = sym2poly(expand(dz_bsf));                            
k = dz_bsf(1);                                            
dz_bsf = dz_bsf/k;
nz_bsf = nz_bsf/k;

%% Cascade Bandpass and Bandstop Filters
% Multiply the transfer functions of both filters
nz_cascade = conv(nz_bpf, nz_bsf);
dz_cascade = conv(dz_bpf, dz_bsf);

% Display Results
disp('Cascaded Filter Coefficients:');
fvtool(nz_cascade, dz_cascade)                                        

% Magnitude Plot (not in log scale)
[H, f] = freqz(nz_cascade, dz_cascade, 1024*1024, f_samp);
plot(f, abs(H))
title('Cascaded Filter Magnitude Response')
xlabel('Frequency (KHz)')
ylabel('Magnitude')
grid on;
