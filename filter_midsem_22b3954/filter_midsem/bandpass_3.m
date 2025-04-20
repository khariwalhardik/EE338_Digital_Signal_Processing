% Butterworth Analog LPF parameters
Wc = 1.08;        % Cut-off frequency
N = 15;           % Order 

% Calculate poles of Butterworth polynomial of degree 15 in the open CLHP
p = zeros(1, N);  % Pre-allocate poles array
for k = 1:N
    theta = pi/2 + (2*k-1)*pi/(2*N);
    p(k) = Wc * (cos(theta) + 1i*sin(theta));
end

% Band Edge Specifications
fp1 = 70;
fs1 = 75;
fs2 = 210;
fp2 = 215;

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
fvtool(nz,dz)                                         % Frequency response using fvtool

% Magnitude Plot (not in log scale)
[H, f] = freqz(nz, dz, 1024*1024, f_samp);
plot(f, abs(H))
title('Magnitude Response')
xlabel('Frequency (KHz)')
ylabel('Magnitude')
grid on;
