%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'TankConnectorMk5_2Dxy1';
path = strcat(folder_path,file_name,'.mat');
load(path)
disp('Data Timestamp:')
disp(scpSettings.timestamp)
% check size of data array - note whether one or two channel
disp('Data Size:')
disp(size(scanData))

%%
scanData = squeeze(scanData(:,:,1,:,:,:));
disp(size(scanData))

t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

%% Average repeats

means = mean(scanData,5);
disp(size(means))

RMS = squeeze(rms(means,3));
disp(size(RMS))

%% Vpk-pk

pkrange = [0,20]; % us - time range to look for peak 
pkrangeidx = ceil(pkrange*scpSettings.SampleFrequency/1e6)+1; % corresponding array index

Vpk = squeeze(max(means(:,:,pkrangeidx(1):pkrangeidx(2),1),[],3));
disp(size(Vpk))

%% Example Raw Waveform

xidx = 1;
yidx = 1;
nthrepeat = 5;

wvfm = squeeze(scanData(xidx,yidx,:,1,nthrepeat));
trig = squeeze(scanData(xidx,yidx,:,2,nthrepeat));
input = squeeze(scanData(xidx,yidx,:,3,nthrepeat));

figure(1)
tiledlayout(3,1)
nexttile
plot(t,wvfm)
xlabel('Time (us)')
ylabel('Amplitude')
nexttile
plot(t,trig) %scaled to match peak of waveform
xlabel('Time (us)')
ylabel('Trigger')
nexttile
plot(t,input)
xlabel('Time (us)')
ylabel('Input')

%% Example Mean Waveform

xidx = 10;
yidx = 10;
nthrepeat = 5;

wvfm = squeeze(means(xidx,yidx,:,1));
trig = squeeze(means(xidx,yidx,:,2));
input = squeeze(means(xidx,yidx,:,3));

figure(2)
tiledlayout(3,1)
nexttile
plot(t,wvfm) 
hold on
xline(pkrange)
hold off
xlabel('Time (us)')
ylabel('Amplitude')
nexttile
plot(t,trig) %scaled to match peak of waveform
xlabel('Time (us)')
ylabel('Trigger')
nexttile
plot(t,input)
xlabel('Time (us)')
ylabel('Input')

%%
RMS_Field = squeeze(RMS(:,:,1));
figure(3)
imagesc(RMS_Field)
colorbar
xlabel('X Index')
ylabel('Y Index')
title('RMS Field')

%%
figure(4)
imagesc(Vpk)
colorbar
xlabel('X Index')
ylabel('Y Index')
title('Vpk-pk')