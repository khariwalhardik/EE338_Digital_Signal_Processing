% Butterworth Analog LPF parameters
Wc = 1.08;        % Cut-off frequency
N = 17;           % Order 

% Calculate poles of Butterworth polynomial of degree N in the open CLHP
p = zeros(1, N);  % Pre-allocate poles array
for k = 1:N
    theta = pi/2 + (2*k-1)*pi/(2*N);
    p(k) = Wc * (cos(theta) + 1i*sin(theta));
end

% Band Edge Specifications
fp1 = 79;
fs1 = 84;
fs2 = 203;
fp2 = 208;

% Transformed Band Edge specs using Bilinear Transformation
f_samp = 630;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

% Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

% Analog LPF Transfer Function
[num, den] = zp2tf([], p, Wc^N);

% Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);      % Analog LPF Transfer Function
analog_bpf(s) = analog_lpf((s*s + W0*W0)/(B*s));      % Bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));            % Bilinear transformation

% Coeffs of analog BPF
[ns, ds] = numden(analog_bpf(s));                     
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          
k = ds(1);    
ds = ds/k;
ns = ns/k;

% Coeffs of discrete BPF
[nz, dz] = numden(discrete_bpf(z));                   
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                            
k = dz(1);                                            
dz = dz/k;
nz = nz/k;

% Frequency Response
%fvtool(nz,dz)                                         % Frequency response using fvtool

% Magnitude Plot (not in log scale)
[H, f] = freqz(nz, dz, 1024*1024, f_samp);
plot(f, abs(H))
title('Magnitude Response')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
grid on;

% Printing Magnitude at specified frequencies
freqs_to_check = [70, 75, 210, 215];
mag_values = abs(interp1(f, H, freqs_to_check));  % Interpolating magnitude values

% Display the results
fprintf('Magnitude at 70 Hz: %.4f\n', mag_values(1));
fprintf('Magnitude at 75 Hz: %.4f\n', mag_values(2));
fprintf('Magnitude at 210 Hz: %.4f\n', mag_values(3));
fprintf('Magnitude at 215 Hz: %.4f\n', mag_values(4));
