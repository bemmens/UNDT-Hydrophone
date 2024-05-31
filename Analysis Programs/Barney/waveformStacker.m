clearvars

folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = strcat('PinCyclePin',string(3));
path = strcat(folder_path,file_name,'.mat');
load(path)

bigness = size(scanData);

size(scanData)
x_index = 13;
y_index = 13;
data = squeeze(scanData(x_index,y_index,:,1));

pinNums = [1,2,3,4,6,7,8,9];
waveforms = zeros(bigness(3),length(pinNums));

folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';

for i = 1:length(pinNums)
    pinNum = pinNums(i);
    file_name = strcat('PinCyclePin',string(pinNum));
    path = strcat(folder_path,file_name,'.mat');
    load(path)
    waveforms(:,i) = squeeze(scanData(x_index,y_index,:,:));;
end
%%
%waveformsflat = reshape(waveforms, [], bigness(3), length(pinNums));

%%
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

figure(1)
for i = 1:length(waveforms(1,:))
    plot(t,waveforms(:,i))
    xlim([25,60])
    hold on
end
xlabel('Time [us]')
ylabel('Voltage [V]')
title('Waveform for Single Pin Activation')
legend(string(pinNums))
hold off

%%
variation = std(waveforms,0,2);
plot(t,variation)
xlim([25,60])
xlabel('Time [us]')
ylabel('Voltage [V]')
title('Standard Deviation of Waveforms Across All Pins')


