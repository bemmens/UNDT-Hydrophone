
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

pkrange = [1,50]; % us - time range to look for peak 
pkrangeidx = pkrange*scpSettings.SampleFrequency/1e6; % corresponding array index
[Vpk,idx] = max(saveData.data(:,pkrangeidx(1):pkrangeidx(2)),[],2);

point = 40;
figure(1)
plot(saveData.timestamps/60,Vpk)
%xline(saveData.timestamps(point)/60)
xlabel('Time [mins]')
ylabel('Voltage [V]')


figure(2)
plot(t,saveData.data(point,:))
xline(t(idx(point)))
xline(t(pkrangeidx(1)))
xline(t(pkrangeidx(2)))