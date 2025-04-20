



%----------------------------------------------------------------------------------------
% Sampling Frequency
f_samp = 630;  

% Tolerances (Maximum Ripple)
delta1 = 0.15;  % Passband ripple
delta2 = 0.15;  % Stopband ripple

% Calculate d1 and d2 using delta1 and delta2
d1 = sqrt((1/(1-delta1)^2) - 1);
d2 = sqrt((1/delta2^2) - 1);

%% Group-I Specifications
% Passband and Stopband Frequencies for Group-I
fp1_g1 = 75; fs1_g1 = 70;
fs2_g1 = 110; fp2_g1 = 105;

% Normalized digital frequencies using Bilinear Transformation
wp1_g1 = tan(fp1_g1 / f_samp * pi);
ws1_g1 = tan(fs1_g1 / f_samp * pi);
ws2_g1 = tan(fs2_g1 / f_samp * pi);
wp2_g1 = tan(fp2_g1 / f_samp * pi);

% Bandwidth and Center Frequency
B_g1 = wp2_g1 - wp1_g1;
w0_g1 = sqrt(wp1_g1 * wp2_g1);

% Calculate Order for Group-I Bandpass Filter
N_g1_1 = ceil(0.5 * log10(d2/d1) / log10(ws1_g1/wp1_g1));
N_g1_2 = ceil(0.5 * log10(d2/d1) / log10(ws2_g1/wp2_g1));
N_bp_g1 = max(N_g1_1, N_g1_2);

% Display Results for Group-I
disp('Group-I Results:');
fprintf('Bandpass Filter Order: N = %d\n', N_bp_g1);

%% Group-II Specifications
% Passband and Stopband Frequencies for Group-II
fp1_g2 = 180; fs1_g2 = 175;
fs2_g2 = 215; fp2_g2 = 210;

% Normalized digital frequencies using Bilinear Transformation
wp1_g2 = tan(fp1_g2 / f_samp * pi);
ws1_g2 = tan(fs1_g2 / f_samp * pi);
ws2_g2 = tan(fs2_g2 / f_samp * pi);
wp2_g2 = tan(fp2_g2 / f_samp * pi);

% Bandwidth and Center Frequency
B_g2 = wp2_g2 - wp1_g2;
w0_g2 = sqrt(wp1_g2 * wp2_g2);

% Calculate Order for Group-II Bandpass Filter
N_g2_1 = ceil(0.5 * log10(d2/d1) / log10(ws1_g2/wp1_g2));
N_g2_2 = ceil(0.5 * log10(d2/d1) / log10(ws2_g2/wp2_g2));
N_bp_g2 = max(N_g2_1, N_g2_2);

% Display Results for Group-II
disp('Group-II Results:');
fprintf('Bandpass Filter Order: N = %d\n', N_bp_g2);

