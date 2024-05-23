%% Fourier Transform
% https://uk.mathworks.com/help/matlab/ref/fft.html?searchHighlight=fft&s_tid=srchtitle_support_results_1_fft

%% Load Data
load test.mat

%% Input Data
signal = squeeze(scanData(4,1,:,1)); % Time series input
Fs = 100e6;            % Sampling frequency (Hz)     

%% Processing
T = 1/Fs;             % Sampling period       
L = length(signal);   % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(signal);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs/L*(0:(L/2));

%% Plots
figure(1)
plot(f,P1,"LineWidth",3) 

xlim([0,1e4])

title("Single-Sided Amplitude Spectrum")
xlabel("Frequency (Hz)")
ylabel("Amplitude")
