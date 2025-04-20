% Sampling Frequency
f_samp = 500e3;

% Band Edge specifications for Bandstop
fs1 = 103e3;
fp1 = 105e3;
fp2 = 180e3;
fs2 = 182e3;

% Kaiser parameters
A = -20*log10(0.15);
if (A < 21)
    beta = 0;
elseif (A < 51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

% Normalized digital frequencies
ws1 = (fs1/f_samp)*pi*2;
wp1 = (fp1/f_samp)*pi*2;
wp2 = (fp2/f_samp)*pi*2;
ws2 = (fs2/f_samp)*pi*2;

% Transition Width
delta_w = min((wp1 - ws1), (ws2 - wp2));

% Minimum order
N_min = ceil((A-7.95) / (2.285*delta_w));  % Transition width
n = N_min + 13;                           % Adding margin

% Ideal bandstop impulse response of length "n"
bs_ideal = ideal_lp(pi, n) - ideal_lp(wp2, n) + ideal_lp(wp1, n);

% Kaiser Window of length "n" with shape parameter beta
kaiser_win = (kaiser(n, beta))';

% FIR Bandstop Filter
FIR_BandStop = bs_ideal .* kaiser_win;

% Frequency Response using fvtool
fvtool(FIR_BandStop);

% Magnitude Response Plot
[H, f] = freqz(FIR_BandStop, 1, 1024, f_samp);
plot(f, abs(H));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Response of FIR Bandstop Filter');
grid;
