% Butterworth Analog LPF parameters
Wc = 1.08;        
N = 17;           

% Calculate poles of Butterworth polynomial of degree N in the open CLHP
p = zeros(1, N);  
for k = 1:N
    theta = pi/2 + (2*k-1)*pi/(2*N);
    p(k) = Wc * (cos(theta) + 1i*sin(theta));
end

% Bandstop Filter Specifications
fp1 = 103.6;
fs1 = 108.6;
fs2 = 176.6;
fp2 = 181.6;
f_samp = 630;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

% Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

% Analog LPF Transfer Function
[num, den] = zp2tf([], p, Wc^N);

% Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);      
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));      
discrete_bsf(z) = analog_bsf((z-1)/(z+1));            

% Coeffs of discrete BSF
[nz_bsf, dz_bsf] = numden(discrete_bsf(z));                   
nz_bsf = sym2poly(expand(nz_bsf));
dz_bsf = sym2poly(expand(dz_bsf));                            
k = dz_bsf(1);                                            
dz_bsf = dz_bsf/k;
nz_bsf = nz_bsf/k;

% Bandpass Filter Specifications
fp1 = 77;
fs1 = 82;
fs2 = 202;
fp2 = 207;

% Transformed Band Edge specs using Bilinear Transformation
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

% Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

% Analog BPF Transfer Function
analog_bpf(s) = analog_lpf((s*s + W0*W0)/(B*s));      
discrete_bpf(z) = analog_bpf((z-1)/(z+1));            

% Coeffs of discrete BPF
[nz_bpf, dz_bpf] = numden(discrete_bpf(z));                   
nz_bpf = sym2poly(expand(nz_bpf));
dz_bpf = sym2poly(expand(dz_bpf));                            
k = dz_bpf(1);                                            
dz_bpf = dz_bpf/k;
nz_bpf = nz_bpf/k;

% Cascading Filters: Convolution of Transfer Functions
nz_cascaded = conv(nz_bsf, nz_bpf);
dz_cascaded = conv(dz_bsf, dz_bpf);

% Frequency Response of Cascaded Filter
[H_cascaded, f] = freqz(nz_cascaded, dz_cascaded, 1024*1024, f_samp);
plot(f, abs(H_cascaded))
title('Magnitude Response of Cascaded Filter')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
grid on;

% Printing Magnitude at specified frequencies for Bandstop
freqs_to_check_bsf = [105, 110, 175, 180];
mag_values_bsf = abs(interp1(f, H_cascaded, freqs_to_check_bsf));  

fprintf('Magnitude at 105 Hz (BSF): %.4f\n', mag_values_bsf(1));
fprintf('Magnitude at 110 Hz (BSF): %.4f\n', mag_values_bsf(2));
fprintf('Magnitude at 175 Hz (BSF): %.4f\n', mag_values_bsf(3));
fprintf('Magnitude at 180 Hz (BSF): %.4f\n', mag_values_bsf(4));

% Printing Magnitude at specified frequencies for Bandpass
freqs_to_check_bpf = [70, 75, 210, 215];
mag_values_bpf = abs(interp1(f, H_cascaded, freqs_to_check_bpf));  

fprintf('Magnitude at 70 Hz (BPF): %.4f\n', mag_values_bpf(1));
fprintf('Magnitude at 75 Hz (BPF): %.4f\n', mag_values_bpf(2));
fprintf('Magnitude at 210 Hz (BPF): %.4f\n', mag_values_bpf(3));
fprintf('Magnitude at 215 Hz (BPF): %.4f\n', mag_values_bpf(4));
