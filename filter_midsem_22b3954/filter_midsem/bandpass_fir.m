% Sampling Frequency
f_samp = 500e3;

% Band Edge specifications
fs1 = 73e3;
fp1 = 75e3;
fp2 = 210e3;
fs2 = 212e3;

% Normalized digital frequencies
ws1 = (fs1/f_samp)*2*pi;
wp1 = (fp1/f_samp)*2*pi;
wp2 = (fp2/f_samp)*2*pi;
ws2 = (fs2/f_samp)*2*pi;

% Kaiser parameters
A = -20*log10(0.15);
if (A < 21)
    beta = 0;
elseif (A < 51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

% Minimum order
N_min = ceil((A-7.95) / (2.285*(wp1-ws1)));  % Transition width wp1-ws1
n = N_min + 13;                             % Adding margin

% Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(wp2, n) - ideal_lp(wp1, n);

% Kaiser Window of length "n" with shape parameter beta
kaiser_win = (kaiser(n, beta))';

% FIR Bandpass Filter
FIR_BandPass = bp_ideal .* kaiser_win;

% Frequency Response using fvtool
fvtool(FIR_BandPass);

% Magnitude Response Plot
[H, f] = freqz(FIR_BandPass, 1, 1024, f_samp);
plot(f, abs(H));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Response of FIR Bandpass Filter');
grid;
