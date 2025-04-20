% Sampling Frequency
f_samp = 1000;

% Bandstop 1: Stops everything except 50 - 100 Hz
fp1_1 = 50; fs1_1 = 45;
fs2_1 = 105; fp2_1 = 100;

% Bandstop 2: Stops everything except 200 - 250 Hz
fp1_2 = 200; fs1_2 = 195;
fs2_2 = 255; fp2_2 = 250;

% Convert to digital frequencies
wp1_1 = tan(fp1_1/f_samp*pi); ws1_1 = tan(fs1_1/f_samp*pi);
ws2_1 = tan(fs2_1/f_samp*pi); wp2_1 = tan(fp2_1/f_samp*pi);

wp1_2 = tan(fp1_2/f_samp*pi); ws1_2 = tan(fs1_2/f_samp*pi);
ws2_2 = tan(fs2_2/f_samp*pi); wp2_2 = tan(fp2_2/f_samp*pi);

% Butterworth Bandstop Design Parameters
Wc = 1.08; N = 5;

% Analog LPF Poles for Butterworth Filter
p = zeros(1, N);  
for k = 1:N
    theta = pi/2 + (2*k-1)*pi/(2*N);
    p(k) = Wc * (cos(theta) + 1i*sin(theta));
end
[num, den] = zp2tf([], p, Wc^N);

syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);

% Bandstop 1 Transformation
W0_1 = sqrt(wp1_1*wp2_1);
B_1 = wp2_1-wp1_1;
analog_bsf1(s) = analog_lpf((B_1*s)/(s*s + W0_1*W0_1)); 
discrete_bsf1(z) = analog_bsf1((z-1)/(z+1));

% Coefficients of Discrete BSF 1
[nz_bs1, dz_bs1] = numden(discrete_bsf1(z));                   
nz_bs1 = sym2poly(expand(nz_bs1));
dz_bs1 = sym2poly(expand(dz_bs1));                            
k_bs1 = dz_bs1(1);                                            
dz_bs1 = dz_bs1/k_bs1;
nz_bs1 = nz_bs1/k_bs1;

% Bandstop 2 Transformation
W0_2 = sqrt(wp1_2*wp2_2);
B_2 = wp2_2-wp1_2;
analog_bsf2(s) = analog_lpf((B_2*s)/(s*s + W0_2*W0_2)); 
discrete_bsf2(z) = analog_bsf2((z-1)/(z+1));

% Coefficients of Discrete BSF 2
[nz_bs2, dz_bs2] = numden(discrete_bsf2(z));                   
nz_bs2 = sym2poly(expand(nz_bs2));
dz_bs2 = sym2poly(expand(dz_bs2));                            
k_bs2 = dz_bs2(1);                                            
dz_bs2 = dz_bs2/k_bs2;
nz_bs2 = nz_bs2/k_bs2;

% Cascade the two Bandstop filters
nz_cascaded = conv(nz_bs1, nz_bs2);
dz_cascaded = conv(dz_bs1, dz_bs2);

% Frequency Response of Bandstop 1
[H_bs1, w_bs1] = freqz(nz_bs1, dz_bs1, 1024, f_samp);
figure;
plot(w_bs1, abs(H_bs1));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Response of Bandstop 1 (50-100 Hz Pass)');
grid on;

% Frequency Response of Bandstop 2
[H_bs2, w_bs2] = freqz(nz_bs2, dz_bs2, 1024, f_samp);
figure;
plot(w_bs2, abs(H_bs2));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Response of Bandstop 2 (200-250 Hz Pass)');
grid on;

% Frequency Response of Cascaded Filter (Multibandpass)
[H_cascaded, w_cascaded] = freqz(nz_cascaded, dz_cascaded, 1024, f_samp);
figure;
plot(w_cascaded, abs(H_cascaded));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Response of Cascaded Multibandpass Filter');
grid on;
