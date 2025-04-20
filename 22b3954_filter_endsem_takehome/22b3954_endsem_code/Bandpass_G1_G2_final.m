clc; clear; close all;
f_samp = 630e3;  % Sampling frequency

% Transition band edges for both bands
% Band 1: 75kHz–105kHz (transition: 70kHz–110kHz)
fs1_1 = 70e3; fp1_1 = 75e3; fp2_1 = 105e3; fs2_1 = 110e3;

% Band 2: 180kHz–210kHz (transition: 175kHz–215kHz)
fs1_2 = 175e3; fp1_2 = 180e3; fp2_2 = 210e3; fs2_2 = 215e3;

% Normalized cutoffs for each band (use midpoints of transitions)
Wc1_1 = ((fs1_1 + fp1_1)/2)*2*pi/f_samp;
Wc2_1 = ((fp2_1 + fs2_1)/2)*2*pi/f_samp;

Wc1_2 = ((fs1_2 + fp1_2)/2)*2*pi/f_samp;
Wc2_2 = ((fp2_2 + fs2_2)/2)*2*pi/f_samp;

fprintf('Band 1 Wc1 = %.3f, Wc2 = %.3f\n', Wc1_1, Wc2_1);
fprintf('Band 2 Wc1 = %.3f, Wc2 = %.3f\n', Wc1_2, Wc2_2);

% Common Kaiser parameters
delta = 0.15;
A = -20*log10(delta);
fprintf('A = %.3f\n', A);

if A < 21
    beta = 0;
elseif A < 51
    beta = 0.5842*(A - 21)^0.4 + 0.07886*(A - 21);
else
    beta = 0.1102*(A - 8.7);
end

% Get the tighter transition width
delta_W = min([Wc2_1 - Wc1_1, Wc2_2 - Wc1_2]);
N_min = ceil((A - 7.95) / (2.285 * delta_W));
n = N_min + 100;
if mod(n, 2) == 0
    n = n + 1;
end
fprintf('Final filter length N = %d\n', n);

% Ideal bandpass impulse responses
bp1_ideal = ideal_lp(Wc2_1, n) - ideal_lp(Wc1_1, n);
bp2_ideal = ideal_lp(Wc2_2, n) - ideal_lp(Wc1_2, n);

% Kaiser window
kaiser_win = (kaiser(n, beta))';

% Windowed filters
FIR_Band1 = bp1_ideal .* kaiser_win;
FIR_Band2 = bp2_ideal .* kaiser_win;

% Final multiband FIR filter (sum of both)
FIR_MultiBand = FIR_Band1 + FIR_Band2;

% Visualize individual bands and combined filter
fvtool(FIR_Band1, FIR_Band2, FIR_MultiBand);
legend('Band 1', 'Band 2', 'Combined MultiBand');

% Frequency response
[H, f] = freqz(FIR_MultiBand, 1, 1024, f_samp);
figure;
plot(f, abs(H)); grid on; hold on;
xline([fs1_1, fp1_1, fp2_1, fs2_1, fs1_2, fp1_2, fp2_2, fs2_2], '--r');
title('Magnitude Response of Combined Multi-Band FIR Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Magnitude Response', 'Band Edges');

% Magnitude at corners
corner_freqs = [fs1_1, fp1_1, fp2_1, fs2_1, fs1_2, fp1_2, fp2_2, fs2_2];
for i = 1:length(corner_freqs)
    [~, idx] = min(abs(f - corner_freqs(i)));
    fprintf('Magnitude at %.1f Hz: %.4f\n', corner_freqs(i), abs(H(idx)));
end
