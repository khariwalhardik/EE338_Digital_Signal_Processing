clc; clear; close all;

f_samp = 630e3;  % Sampling frequency

%% Group-1 Bandpass Filter (75 kHz to 105 kHz)
fs1_1 = 70e3; fp1_1 = 75e3; fp2_1 = 105e3; fs2_1 = 110e3;
Wc1_1 = ((fs1_1 + fp1_1)/2)*2*pi/f_samp;
Wc2_1 = ((fp2_1 + fs2_1)/2)*2*pi/f_samp;

A = -20*log10(0.15);  % Stopband attenuation
if A < 21
    beta = 0;
elseif A < 51
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

wct1 = (fp1_1 - fs1_1)*2*pi/f_samp;
N_min1 = ceil((A - 7.95) / (2.285 * wct1));
n1 = N_min1 + 23;
if mod(n1, 2) == 0
    n1 = n1 + 1;
end

bp_ideal_1 = ideal_lp(Wc2_1, n1) - ideal_lp(Wc1_1, n1);
kaiser_win_1 = kaiser(n1, beta)';
FIR_BandPass_1 = bp_ideal_1 .* kaiser_win_1;

%% Group-2 Bandpass Filter (180 kHz to 210 kHz)
fs1_2 = 175e3; fp1_2 = 180e3; fp2_2 = 210e3; fs2_2 = 215e3;
Wc1_2 = ((fs1_2 + fp1_2)/2)*2*pi/f_samp;
Wc2_2 = ((fp2_2 + fs2_2)/2)*2*pi/f_samp;

wct2 = (fp1_2 - fs1_2)*2*pi/f_samp;
N_min2 = ceil((A - 7.95) / (2.285 * wct2));
n2 = N_min2 + 20;
if mod(n2, 2) == 0
    n2 = n2 + 1;
end

bp_ideal_2 = ideal_lp(Wc2_2, n2) - ideal_lp(Wc1_2, n2);
kaiser_win_2 = kaiser(n2, beta)';
FIR_BandPass_2 = bp_ideal_2 .* kaiser_win_2;

%% Combine using Parallel Structure (Zero pad to match lengths)
len_max = max(length(FIR_BandPass_1), length(FIR_BandPass_2));
FIR_BandPass_1 = [FIR_BandPass_1, zeros(1, len_max - length(FIR_BandPass_1))];
FIR_BandPass_2 = [FIR_BandPass_2, zeros(1, len_max - length(FIR_BandPass_2))];

FIR_BandPass = FIR_BandPass_1 + FIR_BandPass_2;





disp('FIR Filter Coefficients(Final-filter):');
disp(FIR_BandPass);



%% ------------------ FIGURE 1 ------------------
% Frequency response (both linear and log scale)
figure(1);
freqz(FIR_BandPass, 1, 1024, f_samp);
title('Figure 1: Frequency Response (Discrete & Log Scale)');

%% ------------------ FIGURE 2 ------------------
% Magnitude response (unnormalized)
[H, f] = freqz(FIR_BandPass, 1, 1024, f_samp);
figure(2); clf;
plot(f, abs(H), 'LineWidth', 1.5); grid on; hold on;

% Define band edge frequencies
corner_frequencies = [fs1_1, fp1_1, fp2_1, fs2_1,fs1_2,fp1_2,fp2_2,fs2_2];
xline(corner_frequencies, '--r');

% Add frequency and magnitude values as data-tip style labels
for i = 1:length(corner_frequencies)
    % Find nearest frequency index
    [~, idx] = min(abs(f - corner_frequencies(i)));
    freq_val = f(idx);
    mag_val = abs(H(idx));
    
    % Plot black dot at that location
    plot(freq_val, mag_val, 'ko', 'MarkerFaceColor', 'k');

    % Add text box near the point
    text(freq_val, mag_val, ...
        sprintf('X %.0f\nY %.4f', freq_val, mag_val), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 9, ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'Margin', 4);
end

title('Figure 2: Magnitude Response (Unnormalized)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Magnitude Response', 'Band Edges');


% Print magnitudes at band edges
corner_frequencies = [fs1_1, fp1_1, fp2_1, fs2_1,fs1_2,fp1_2,fp2_2,fs2_2];
for i = 1:length(corner_frequencies)
    [~, idx] = min(abs(f - corner_frequencies(i)));
    fprintf('Magnitude at %.1f Hz: %.4f\n', corner_frequencies(i), abs(H(idx)));
end

%% ------------------ FIGURE 3 ------------------
% Phase response
figure(3);
plot(f, unwrap(angle(H)), 'LineWidth', 1.5); grid on;
title('Figure 3: Phase Response (Unnormalized)');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');

%% ------------------ FIGURE 4 ------------------
% Impulse response
figure(4);
stem(FIR_BandPass, 'filled');
title('Figure 4: Impulse Response of Bandpass Filter');
xlabel('n');
ylabel('h[n]');

%% ------------------ FIGURE 5 ------------------
% Coefficients bar graph
figure(5);
bar(FIR_BandPass);
title('Figure 5: FIR Filter Coefficients');
xlabel('Index');
ylabel('Coefficient Value');

%% ------------------ FIGURE 6 ------------------
% Zoomed-in coefficients (optional for better visibility)
figure(6);
plot(FIR_BandPass, 'LineWidth', 1.2);
title('Figure 6: Zoomed Coefficients of Bandpass Filter');
xlabel('Index');
ylabel('Coefficient Value');
grid on;