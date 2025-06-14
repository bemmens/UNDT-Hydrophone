%% Fourier Transform
% https://uk.mathworks.com/help/matlab/ref/fft.html?searchHighlight=fft&s_tid=srchtitle_support_results_1_fft

%% Load Data
% Consider clearing workspace to relieve RAM
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = '419kHz_11';                                                                                         %CHECK
path = strcat(folder_path,file_name,'.mat');
load(path)

%% Input Data
size(scanData)
x_index = 5;
y_index = 5;
signal = squeeze(scanData(x_index,y_index,:,1)); % Time series input
Fs = scpSettings.SampleFrequency;            % Sampling frequency (Hz)                                                                      %CHECK

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
plot(f/1e6,P1,"LineWidth",1) 
xlim([0.5,2])                                                                                                       %CHECK
xline(1.816)
title("Single-Sided Amplitude Spectrum")
title('FFT of Pin 9 Waveform')
xlabel("Frequency (MHz)")
ylabel("Amplitude")
yticks([])