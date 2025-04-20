syms s w0 B;
% Low-pass filter denominators
D1 = s^5 - 1.724*s^4 + 2.944*s^3 - 2.611*s^2 + 1.655*s - 0.4899;
D2 = s^5 - 1.357*s^4 + 1.824*s^3 - 1.273*s^2 + 0.6349*s - 0.1479;

% Frequency transformation: s -> (s^2 + w0^2) / (B s)
s_sub = (s^2 + w0^2) / (B * s);

% Transform the denominators
D1_bp = subs(D1, s, s_sub);
D2_bp = subs(D2, s, s_sub);

% Expand the expressions for better readability
D1_bp = expand(D1_bp);
D2_bp = expand(D2_bp);

% Low-pass to bandpass filter transfer function
K1 = 1.469;
K2 = 0.4437;
H_BPF1 = K1 / D1_bp;
H_BPF2 = K2 / D2_bp;

% Display results
disp('Bandpass Filter Transfer Function (Group 1):');
pretty(H_BPF1)
disp('Bandpass Filter Transfer Function (Group 2):');
pretty(H_BPF2)
