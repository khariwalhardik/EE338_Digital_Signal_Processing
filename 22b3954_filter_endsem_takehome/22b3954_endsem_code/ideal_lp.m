function hd = ideal_lp(wc, M)
% ideal_lp generates the impulse response of an ideal lowpass filter
% wc: cutoff frequency in radians (0 to pi)
% M: filter length (should be odd for symmetry)

alpha = (M - 1) / 2;
n = 0:M-1;
m = n - alpha;

% Avoid division by zero at the center (m = 0)
hd = zeros(1, M);
hd(m == 0) = wc / pi;                % sinc(0) = 1, so value at center = wc/pi
hd(m ~= 0) = sin(wc * m(m ~= 0)) ./ (pi * m(m ~= 0));
end
