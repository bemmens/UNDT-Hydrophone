% HS5_DualInstrument_AcquistionScript_v2.m
% written by F.Cegla March2017
% based on Tiepie LibTiepie examples
% script sents and receives a Hanning windowed Toneburst with user
% controllable settings
% ========================================================================
tic
%% Acquistions settings
% (the only variables that have to be set by the user
nRecLength = 96000;
nSamplingFrequency = 50e6;
nAmplitude = 10;
nFrequency = 2e6;
nCycles = 4;
nAverages = 10;
%pathname = 'D:\LaptopData\Research\Projects\EMATdev';
%filename ='OldBoxPowersupply28Oct.mat';


%%

% Open LibTiePie and display library info if not yet opened:
LibTiePieNeeded

% Search for devices:
LibTiePie.DeviceList.update();

% Try to open a generator with arbitrary support:
clear gen1;
clear gen2;

for k = 0 : LibTiePie.DeviceList.Count - 1
    item = LibTiePie.DeviceList.getItemByIndex(k);
%    if k==1
        if item.canOpen(DEVICETYPE.GENERATOR)
            if exist('gen1')
              gen2 = item.openGenerator();
              if ~ismember(GM.BURST_COUNT, gen2.ModesNative)
                clear gen2;
              end
            else
              gen1 = item.openGenerator();
              if ~ismember(GM.BURST_COUNT, gen1.ModesNative)
                clear gen1;
              end
            end  
        end
%    end
%    if k==2
%        if item.canOpen(DEVICETYPE.GENERATOR)
%            gen2 = item.openGenerator();
%            if ~ismember(GM.BURST_COUNT, gen2.ModesNative)
%              clear gen2;
%            end
%        end
%    end
    
end
clear item

%%
% configure Generators

% Generator 
if exist('gen1', 'var')
    % Generator1
    % Set signal type:
    gen1.SignalType = ST.ARBITRARY;
    % Set frequency:
    gen1.FrequencyMode = FM.SAMPLEFREQUENCY;
    gen1.Frequency = nSamplingFrequency; % in Hz
    % Set amplitude:
    gen1.Amplitude = nAmplitude; % in Volts
    % Set offset:
    gen1.Offset = 0; % 0 V
    % Set mode:
    gen1.Mode = GM.BURST_COUNT;
    % Set burst count:
    gen1.BurstCount = 1; % 1 period
    gen1.AmplitudeAutoRanging=1;
    % Enable output:
    gen1.OutputOn = true;
    % Print generator info:
    display(gen1);       
    %%
    % Generator2
    % Set signal type:
    gen2.SignalType = ST.ARBITRARY;
    % Set frequency:
    gen2.FrequencyMode = FM.SAMPLEFREQUENCY;
    gen2.Frequency = nSamplingFrequency; % in Hz
    % Set amplitude:
    gen2.Amplitude = nAmplitude; % in Volts
    % Set offset:
    gen2.Offset = 0; % 0 V
    % Set mode:
    gen2.Mode = GM.BURST_COUNT;
    % Set burst count:
    gen2.BurstCount = 1; % 1 period
    gen2.AmplitudeAutoRanging=1;
    % Enable output:
    gen2.OutputOn = true;
    % Print generator info:
    display(gen2);       
    %%    
    
% create Oscilloscope
 clear scp;
for k = 0 : LibTiePie.DeviceList.Count - 1
    item = LibTiePie.DeviceList.getItemByIndex(k);
    if item.canOpen(DEVICETYPE.OSCILLOSCOPE)
        scp = item.openOscilloscope();
        if ismember(MM.BLOCK, scp.MeasureModes)
            break;
        else
            clear scp;
        end
    end
end
clear item

if exist('scp', 'var')
    % Set measure mode:
    scp.MeasureMode = MM.BLOCK;
    % Set sample frequency:
    scp.SampleFrequency = nSamplingFrequency; % Hz
    % Set record length:
    scp.RecordLength = nRecLength; % nRecLength Samples
    % Set pre sample ratio:
    scp.PreSampleRatio = 0; % 0 %
    % For all channels:
    for ch = scp.Channels
        % Enable channel to measure it:
        ch.Enabled = true;
        % Set range:
        ch.Range = nAmplitude; % Volts
        % Set coupling:
        ch.Coupling = CK.DCV; % DC Volt
        % Release reference:
        clear ch;
    end

    % Set trigger timeout:
    scp.TriggerTimeOut = 0.1;%10e-3; % 100 ms
    % Disable all channel trigger sources:
    for ch = scp.Channels
        ch.Trigger.Enabled = false;
        clear ch;
    end
    % Setup channel trigger:
    chTr = scp.Channels(1).Trigger; % Ch 1
    % Enable trigger source:
    chTr.Enabled = true;
    % Kind:
    chTr.Kind = TK.RISINGEDGE; % Rising edge
    % Level:
    chTr.Levels(1) = 0.1; % 50 %
    % Hysteresis:
    chTr.Hystereses(1) = 0.05; % 5 %
    % Release reference:
    clear chTr;

  %%  
    % Locate trigger input:
    %triggerInput = scp.getTriggerInputById(TIID.GENERATOR_STOP); % or TIID.GENERATOR_START or TIID.GENERATOR_STOP
    %triggerInput = scp.getTriggerInputById(35651585) %genstop 32709
    triggerInput = scp.getTriggerInputById(35651584) %genstart 32709
    %triggerInput = scp.getTriggerInputById(18874368) %genstart 32392
    
    %triggerInput = scp.getTriggerInputById(TIID.GENERATOR_START); % or TIID.GENERATOR_START or TIID.GENERATOR_STOP

    if triggerInput == false
        clear triggerInput;
        clear scp;
        clear gen1;
        error('Unknown trigger input!');
    end
    % Enable trigger input:
    triggerInput.Enabled = true;
    % Release reference to trigger input:
    clear triggerInput;
%%
    % Print oscilloscope info:
    display(scp);

    % Start measurement:
    %scp.start();

    % Wait for measurement to complete:
    %while ~scp.IsDataReady
    %    pause(10e-3) % 10 ms delay, to save CPU time.
    %end

    else
    error('No oscilloscope available with block measurement support!');
end   
 

    
%% create signal to load     
    buffer=zeros(1,scp.RecordLength);
    f=nFrequency;
    t=[0:1/nSamplingFrequency:(scp.RecordLength-1)/nSamplingFrequency];
    TotalpulseLength=floor(nCycles*nSamplingFrequency/nFrequency);
    buffer=zeros(1,TotalpulseLength);
    y=hann(TotalpulseLength)'.*(sin(2.*pi.*f.*t(1:TotalpulseLength)));
   
    
    for count=1:nAverages
        % load signal into generator
        gen1.setData(y);
        %gen2.setData(-y);
        % Start signal burst:
        scp.start();
        gen1.start();
        %gen2.start();
        % Wait for burst to complete:
        while gen1.IsBurstActive
            pause(10e-3); % 10 ms delay, to save CPU time.
            display('gen not ready')
        end
        while ~scp.IsDataReady
            pause(10e-3) % 10 ms delay, to save CPU time.
            display('osc not ready')
        end
        
        % check if measurement was triggered
        Trig=scp.IsTriggered;
        % or if trigger time out was invoked
        TrigTO=scp.IsTimeOutTriggered;
        
        
        %% Get data:
        arData = scp.getData();
        
        % Get all channel data value ranges (which are compensated for probe gain/offset):
        clear darRangeMin;
        clear darRangeMax;
        for i = 1 : length(scp.Channels)
            [darRangeMin(i), darRangeMax(i)] = scp.Channels(i).getDataValueRange();
        end
        if count==1
            Data=arData;
        else
            Data=Data+arData;
        end
    end
    test=Data;
    Data=Data./nAverages;
    %%
   
    % Disable output:
    gen1.OutputOn = false; 
    % Close generator:
    clear gen1;
else
    error('No generator available with burst support!');
end

%remove DC offset
result=double(Data(:,1));
result=xcorr(result,y);
result=result(end-length(Data(:,1))+1:end);
%FR=fft(result);
%FR(1)=0;
%result=real(ifft(FR));
offset=0.0002225-0.0002;
offset=0;
    % Plot Xcorr result:
%     figure;
%     normval=max(xcorr(y,y));
%     
%     plot(((1:length(Data(:,1))) / nSamplingFrequency-offset), result./normval);
%     %plot(((1:length(Data(:,1))) / nSamplingFrequency-offset), result./max(result));
%     %plot(result)
%     %plot((1:length(Data(:,1))) / scp.SampleFrequency, Data./max(Data(:,1)));
%     %axis([0 (scp.RecordLength / scp.SampleFrequency) min(darRangeMin) max(darRangeMax)]);
%     xlabel('Time [s]');
%     ylabel('Amplitude [V]');
%     %ylim([-0.1 0.1]);
    %xlim([0 0.0001])
    
    %plot raw data
    figure
    plot(((1:length(Data(:,1))) / nSamplingFrequency-offset), Data(:,:));
     xlabel('Time [s]');
    ylabel('Amplitude [V]');
    %ylim([-0.1 0.1]);
    %xlim([0 0.0001])
    % Close oscilloscope:
    clear scp;

%% Save data
% fred=datetime;
% fv=datevec(fred);
% FDS=strcat(num2str(fv(1)),'-',num2str(fv(2)),'-',num2str(fv(3)),'-',...
%     num2str(fv(4)),'-',num2str(fv(5)),'-',num2str(floor(fv(6))),'-');
% save(strcat(FDS,filename));
%cd(pathname);
%save(filename);
toc