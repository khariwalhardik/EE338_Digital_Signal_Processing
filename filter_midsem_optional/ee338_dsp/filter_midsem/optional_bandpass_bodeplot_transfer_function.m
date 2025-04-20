clear;
clc;
close all;

%% Group 1: Bandpass Transfer Function
B1 = 0.1846;  
w0_1 = 1.08;  

% Lowpass Transfer Function for Group 1 (Prototype)
num_lp1 = 1.469;  % Given numerator (adjust if necessary)
den_lp1 = [1, -1.357, 1.824, -1.273, 0.6349, -0.1479];  % Given denominator

% Convert Lowpass TF to Bandpass
[num_bp1, den_bp1] = lp2bp(num_lp1, den_lp1, w0_1, B1);

% Create Transfer Function
H_bp1 = tf(num_bp1, den_bp1);

% Display Transfer Function
disp('Corrected Bandpass Transfer Function (Group 1):');
disp(H_bp1);

% Plot Frequency Response
figure;
bode(H_bp1);
title('Corrected Bode Plot of Bandpass Filter (Group 1)');
grid on;

%% Group 2: Bandpass Transfer Function
B2 = 0.3278;  
w0_2 = 2.50;  

% Lowpass Transfer Function for Group 2 (Prototype)
num_lp2 = 0.4437;  % Given numerator (adjust if necessary)
den_lp2 = [1, -1.357, 1.824, -1.273, 0.6349, -0.1479];  % Given denominator

% Convert Lowpass TF to Bandpass
[num_bp2, den_bp2] = lp2bp(num_lp2, den_lp2, w0_2, B2);

% Create Transfer Function
H_bp2 = tf(num_bp2, den_bp2);

% Display Transfer Function
disp('Corrected Bandpass Transfer Function (Group 2):');
disp(H_bp2);

% Plot Frequency Response
figure;
bode(H_bp2);
title('Corrected Bode Plot of Bandpass Filter (Group 2)');
grid on;
