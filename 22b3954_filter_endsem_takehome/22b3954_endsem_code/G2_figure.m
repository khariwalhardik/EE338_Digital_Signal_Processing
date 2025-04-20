clc; clear; close all;

f_samp = 630e3;  % Sampling frequency

% Band Edge specifications for Bandpass
fs1 = 175e3;
fp1 = 180e3;
fp2 = 210e3;
fs2 = 215e3;

% Digital (rad/sec) cutoff frequencies
Wc1 = ((fs1 + fp1)/2)*2*pi/f_samp;
Wc2 = ((fp2 + fs2)/2)*2*pi/f_samp;

fprintf('Wc1 = %.3f \n',Wc1);
fprintf('Wc2 = %.3f \n',Wc2);

% Kaiser window parameters
A = -20*log10(0.15);  % Stopband attenuation
fprintf('A = %.3f \n',A);

if A < 21
    beta = 0;
elseif A < 51
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

wct=(fp1-fs1)*2*pi/f_samp;
N_min = ceil((A-7.95) / (2.285*wct));  % Estimate filter order
n = N_min + 20;                             % Increase length for better response
if mod(n, 2) == 0
    n = n + 1;  % Make it odd
end
fprintf('N = %d \n',n);

% Ideal bandpass impulse response
bp_ideal = ideal_lp(Wc2, n) - ideal_lp(Wc1, n);

% Apply Kaiser window
kaiser_win = kaiser(n, beta)';
FIR_BandPass = bp_ideal .* kaiser_win;

disp('FIR Filter Coefficients(Group-2):');
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
corner_frequencies = [fs1, fp1, fp2, fs2];
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
corner_frequencies = [fs1, fp1, fp2, fs2];
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
