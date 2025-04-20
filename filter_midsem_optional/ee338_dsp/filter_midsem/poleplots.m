% Filter specifications
fs = 250e3;  % Sampling frequency (Hz)
f1 = 75e3;   % Lower passband frequency (Hz)
f2 = 105e3;  % Upper passband frequency (Hz)
tw = 5e3;    % Transition width (Hz)

% Calculate stopband frequencies
fstop1 = f1 - tw/2;
fstop2 = f2 + tw/2;

% Normalize frequencies
wp = [f1 f2] / (fs/2);
ws = [fstop1 fstop2] / (fs/2);

% Ripple specifications
delta1 = 0.15;  % Passband ripple
delta2 = 0.15;  % Stopband attenuation

% Design filter
[n, wn] = buttord(wp, ws, delta1, delta2);
[b, a] = butter(n, wn, 'bandpass');

% Plot frequency response
[h, w] = freqz(b, a, 1024, fs);
plot(w, 20*log10(abs(h)))
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Bandpass Filter Frequency Response')
grid on

% Plot pole-zero diagram
figure
zplane(b, a)
title('Pole-Zero Plot of Bandpass Filter')
