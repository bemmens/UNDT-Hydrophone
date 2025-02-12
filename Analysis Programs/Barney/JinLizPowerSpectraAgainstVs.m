data2load = ["419kHz_2_14_30V","419kHz_2_14_25V","419kHz_2_14_20V","419kHz_2_14_15V","419kHz_2_14_9V"];
dataStructs = struct();

dataStructs.V30 = load(data2load(1));
dataStructs.V25 = load(data2load(2));
dataStructs.V20 = load(data2load(3));
dataStructs.V15 = load(data2load(4));
dataStructs.V9 = load(data2load(5));

Vs = [9,15,20,25,30];
wvfms = [dataStructs.V9.saveData.wvfm,dataStructs.V15.saveData.wvfm,dataStructs.V20.saveData.wvfm,dataStructs.V25.saveData.wvfm,dataStructs.V30.saveData.wvfm];
fs = [dataStructs.V9.saveData.bandpass(3),dataStructs.V15.saveData.bandpass(3),dataStructs.V20.saveData.bandpass(3),dataStructs.V25.saveData.bandpass(3),dataStructs.V30.saveData.bandpass(3)];

Ps = [];
figure;
hold on;
for i = 1:5
    [P, F] = pspectrum(wvfms(:,i), fs(i));
    plot(F, 10*log10(P), 'DisplayName', ['V = ' num2str(Vs(i)) 'V']);
    Ps(i) = P(345);
end
hold off;
xlim([0.3,0.5]*1e6)
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectra');
legend show;
grid on;

figure
plot(Vs,pow2db(Ps))
xlabel('Voltage (V)');
ylabel('Power (dB)');
title('Power vs Voltage')

%%
Fs = scpSettings.SampleFrequency;
low = 0.4*1e6;
high = 0.44*1e6;
b_1 = bandpass(wvfms(:,1)*1e3,[low high],Fs);
b_end = bandpass(wvfms(:,end)*1e3,[low high],Fs);

tiledlayout(2,1)
nexttile
plot(t,wvfms(:,1)*1e3-mean(wvfms(:,1)*1e3))
hold on
plot(t,b_1)
yline(rms(b_1),'-r','bandpassed RMS')
yline(rms(wvfms(:,1)*1e3-mean(wvfms(:,1)*1e3)),'-b','raw RMS')
hold off
xlim([0,50])
xlabel('Time [us]')
ylabel('Volts [mV]')
title('420kHz Signal Off - Surface Not Perturbed')
%legend('raw data','bandpassed','bandpassed RMS',c)
ylim([-30,30])

nexttile
plot(t,wvfms(:,end)*1e3-mean(wvfms(:,end)*1e3))
hold on
plot(t,b_end)
yline(rms(b_end),'-r','bandpassed RMS')
yline(rms(wvfms(:,end)*1e3-mean(wvfms(:,end)*1e3)),'-b','raw RMS')
hold off
xlim([0,50])
xlabel('Time [us]')
ylabel('Volts [mV]')
title('420kHz Signal On - Surface Perturbed')
ylim([-30,30])


%%

bandpass(wvfms(1:1000,end)*1e3,[low high],Fs)