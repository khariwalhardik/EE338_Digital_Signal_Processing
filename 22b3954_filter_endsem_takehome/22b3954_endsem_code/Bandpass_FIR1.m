clc; clear; close all;
f_samp = 630e3;

% Band Edge specifications
fs1 = 70e3;
fp1 = 75e3;
fp2 = 105e3;
fs2 = 110e3;

Wc1 = ((fs1 + fp1)/2)*2*pi/f_samp;
Wc2 = ((fp2 + fs2)/2)*2*pi/f_samp;
fprintf('Wc1 = %.3f \n',Wc1);
fprintf('Wc2 = %.3f \n',Wc2);


% Kaiser parameters
A = -20*log10(0.15);
fprintf('A = %.3f \n',A);
if(A < 21)
    beta = 0;
elseif(A < 51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-7.95) / (2.285*(Wc2-Wc1)));  % Empirical formula for N_min

% Window length for Kaiser Window
n = N_min + 100;
if mod(n, 2) == 0

    n = n + 1;
end
fprintf('N = %d \n',n);
% Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(Wc2,n) - ideal_lp(Wc1,n);

% Kaiser Window of length "n" with shape parameter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass); % Frequency response


%magnitude response
[H,f] = freqz(FIR_BandPass,1,1024, f_samp);
plot(f, abs(H)); grid on; hold on;
xline([fs1, fp1, fp2, fs2], '--r');
title('Magnitude Response with Band Edges');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Magnitude Response', 'Band Edges');

% Find magnitudes at corner frequencies
corner_frequencies = [fs1, fp1, fp2, fs2];
for i = 1:length(corner_frequencies)
    [~, idx] = min(abs(f - corner_frequencies(i))); % Find closest frequency index
    fprintf('Magnitude at %.1f Hz: %.4f\n', corner_frequencies(i), abs(H(idx)));
end
