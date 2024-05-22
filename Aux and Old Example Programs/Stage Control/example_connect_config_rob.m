function frq2 = example_connect_config_rob(device_id, freq, amplitudeOut, no_harms)
% EXAMPLE_CONNECT_CONFIG Connect to and configure a Zurich Instruments device
%
% USAGE R = EXAMPLE_CONNECT_CONFIG(DEVICE_ID)
%  
% Connect to the Zurich Instruments instrument specified by DEVICE_ID, obtain
% a single demodulator sample and calculate its RMS amplitude R. DEVICE_ID
% should be a string, e.g., 'dev2006' or 'uhf-dev2006'.
%
% NOTE Additional configuration: Connect signal output 1 to signal input 1
% with a BNC cable.
%  
% NOTE This is intended to be a simple example demonstrating how to connect
% to a Zurich Instruments device from ziPython. In most cases, data acquistion
% should use either ziDAQServer's poll() method or an instance of the
% ziDAQRecorder class.
%
% NOTE Please ensure that the ziDAQ folders 'Driver' and 'Utils' are in your
% Matlab path. To do this (temporarily) for one Matlab session please navigate
% to the ziDAQ base folder containing the 'Driver', 'Examples' and 'Utils'
% subfolders and run the Matlab function ziAddPath().  
% >>> ziAddPath;
%
% Use either of the commands:
% >>> help ziDAQ
% >>> doc ziDAQ 
% in the Matlab command window to obtain help on all available ziDAQ commands.
%
% See also EXAMPLE_CONNECT, EXAMPLE_POLL.
%
% Copyright 2008-2015 Zurich Instruments AG

addpath(genpath('C:\Program Files\Zurich Instruments\LabOne\API\MATLAB2012'))

clear ziDAQ;

if ~exist('device_id', 'var')
    error(['No value for device_id specified. The first argument to the ' ...
           'example should be the device ID on which to run the example, ' ...
           'e.g. ''dev2006'' or ''uhf-dev2006''.'])
end

% Check the ziDAQ MEX (DLL) and Utility functions can be found in Matlab's path.
if ~(exist('ziDAQ') == 3) && ~(exist('ziDevices', 'file') == 2)
    fprintf('Failed to either find the ziDAQ mex file or ziDevices() utility.\n')
    fprintf('Please configure your path using the ziDAQ function ziAddPath().\n')
    fprintf('This can be found in the API subfolder of your LabOne installation.\n');
    fprintf('On Windows this is typically:\n');
    fprintf('C:\\Program Files\\Zurich Instruments\\LabOne\\API\\MATLAB2012\\\n');
    return
end

% Determine the device identifier from it's ID.
device = lower(ziDAQ('discoveryFind', device_id));

% Get the device's default connectivity properties.
props = ziDAQ('discoveryGet', device);

% The maximum API level supported by this example.
apilevel_example = 5;
% The maximum API level supported by the device class, e.g., MF.
apilevel_device = props.apilevel;
% Ensure we run the example using a supported API level.
apilevel = min(apilevel_device, apilevel_example);
% See the LabOne Programming Manual for an explanation of API levels.

% Create a connection to a Zurich Instruments Data Server (a API session)
% using the device's default connectivity properties.
ziDAQ('connect', props.serveraddress, props.serverport, apilevel);

% Check that the device is visible to the Data Server.
if ~ismember(device, ziDevices())
    message = ['The specified device `', device, '` is not visible to the Data ', ...
               'Server. Please ensure the device is connected by using the LabOne ', ...
               'User Interface or ziControl (HF2 Instruments).'];
    error(message);
end

% Get the device type and its options (in order to set correct device-specific
% configuration).
devtype = ziDAQ('getByte', ['/' device '/features/devtype']);
options = ziDAQ('getByte', ['/' device '/features/options']);

fprintf('Will run the example on `%s`, an `%s` with options `%s`.\n', device, ...
        devtype, regexprep(options, '\n', '|'));

sigout_amplitude = amplitudeOut;
% r1 = run_example(device,devtype,options,sigout_amplitude, freq);
[r2, frq2] = ziSetup(device, devtype, options, sigout_amplitude, freq, no_harms);
% test_example(r,sigout_amplitude);
  
end

function [r2, frq2] = ziSetup(device, devtype, options, sigout_amplitude, freq, no_harms)

fun_frqs = length(freq);
no_frqs = fun_frqs.*no_harms;
if fun_frqs > 2
    message = ['The specified number of frequency fundamentals `', num2str(fun_frqs), '` exceeds number possible.', ...
               'Please ensure the number of fundamental frequencies does not exceed 2.'];
    error(message);
elseif fun_frqs == 2
    
end
if length(sigout_amplitude) == no_frqs

else
    message = ['The specified number of frequency fundamentals `', num2str(fun_frqs), '` exceeds number possible.', ...
               'Please ensure the number of fundamental frequencies does not exceed 2.'];
    error(message);
end

% RUN_EXAMPLE run the example on the specified device
demod_c = '0'; % demod channel
out_c = '0'; % signal output channel
in_c = '0'; % signal input channel
osc_c = '0'; % oscillator

% define the output mixer channel based on the device type and its options
if strfind(devtype,'UHF') & isempty(strfind(options,'MF'))
    out_mixer_c = '3';
elseif strfind(devtype,'HF2LI') & isempty(strfind(options,'MF'))
    out_mixer_c = '6';
elseif strfind(devtype,'MFLI') & isempty(strfind(options,'MD'))
    out_mixer_c = '1';
else
% instrument with multi-frequency option
    out_mixer_c = '0';
end

time_constant = 1e-5; % [s]
demod_rate = 2e3;
% create a base configuration: disable all outputs, demods and scopes
ziDAQ('setDouble',['/' device '/demods/*/rate'], 0.0);
ziDAQ('setInt',['/' device '/demods/*/trigger'], 0);
ziDAQ('setInt',['/' device '/sigouts/*/enables/*'], 0);
if strfind(devtype,'UHF')
% if the device is a UHF additionally disable all demodulators
    ziDAQ('setInt',['/' device '/demods/*/enable'], 0);
    ziDAQ('setInt',['/' device '/scopes/*/enable'], 0);
elseif strfind(devtype,'HF2')
    ziDAQ('setInt',['/' device '/scopes/*/trigchannel'],-1)
end

% configure the device ready for this experiment
ziDAQ('setInt',['/' device '/sigins/' in_c '/imp50'], 0);
ziDAQ('setInt',['/' device '/sigins/' in_c '/ac'], 1);
ziDAQ('setDouble',['/' device '/sigins/' in_c '/range'], 1.2);
ziDAQ('setInt',['/' device '/sigouts/' out_c '/on'], 1);

% Configure the number of outmixer channels that are on
ziDAQ('setDouble',['/' device '/sigouts/' out_c '/range'], 1);
ziDAQ('setDouble',['/' device '/sigouts/' out_c '/amplitudes/*'], 0);
for kk = 0:1:no_frqs-1;
    % RUN_EXAMPLE run the example on the specified device
    out_mixer_c = num2str(kk); % out mixer channel
    ziDAQ('setDouble',['/' device '/sigouts/' out_c '/amplitudes/' out_mixer_c], sigout_amplitude(kk+1)); % outmixer amplitude
    ziDAQ('setDouble',['/' device '/sigouts/' out_c '/enables/' out_mixer_c], 1);
end

if strfind(devtype,'HF2')
    ziDAQ('setInt',['/' device '/sigins/' in_c '/diff'], 0);
    ziDAQ('setInt',['/' device '/sigouts/' out_c '/add'], 0);
end
ziDAQ('setDouble',['/' device '/demods/*/phaseshift'], 0);
ziDAQ('setInt',['/' device '/demods/*/order'], 8);

% Turn on Demods
demod_c = num2str(0);
for harm = 1:1:no_harms;
    for jj = 1:1:length(freq)
        osc_c = num2str(jj-1);
        ziDAQ('setDouble',['/' device '/demods/' demod_c '/rate'], demod_rate);
        ziDAQ('setInt',['/' device '/demods/' demod_c '/harmonic'], harm);
        
        if strfind(devtype,'UHF')
            ziDAQ('setInt',['/' device '/demods/' demod_c '/enable'], 1);
        end
        if strfind(options,'MF')
            % HF2IS and HF2LI multi-frequency option do not support the node oscselect.
            ziDAQ('setInt',['/' device '/demods/' demod_c '/oscselect'], str2double(osc_c));
            ziDAQ('setInt',['/' device '/demods/' demod_c '/adcselect'], str2double(in_c));
        end
        demod_c = num2str(str2double(demod_c)+1);
    end
    ziDAQ('setDouble',['/' device '/demods/*/timeconstant'], time_constant);
    ziDAQ('setDouble',['/' device '/oscs/' osc_c '/freq'], freq(jj)); % [Hz]
end

% Unsubscribe from any streaming data
ziDAQ('unsubscribe','*');

% Pause to get a settled demodulator filter
pause(10*time_constant);

% Perform a global synchronisation between the device and the data server:
% Ensure that the settings have taken effect on the device before issuing the
% ``getSample`` command. Note: ``sync`` must be issued after waiting for the
% demodulator filter to settle above.
ziDAQ('sync');

% obtain a single demodulator sample
for ff =  1:1:no_frqs
    demod_c = num2str(ff-1);
    sample = ziDAQ('getSample',['/' device '/demods/' demod_c '/sample']);
    r2(ff) = sqrt( sample.x^2 + sample.y^2 );
    frq2(ff) = sample.frequency;
    fprintf('Measured RMS amplitude is %fV.\n',r2(ff))
end

end

function r = run_example(device, devtype, options, sigout_amplitude, freq)
% RUN_EXAMPLE run the example on the specified device
demod_c = '0'; % demod channel
out_c = '0'; % signal output channel
in_c = '0'; % signal input channel
osc_c = '0'; % oscillator

% define the output mixer channel based on the device type and its options
if strfind(devtype,'UHF') & isempty(strfind(options,'MF'))
    out_mixer_c = '3';
elseif strfind(devtype,'HF2LI') & isempty(strfind(options,'MF'))
    out_mixer_c = '6';
elseif strfind(devtype,'MFLI') & isempty(strfind(options,'MD'))
    out_mixer_c = '1';
else
% instrument with multi-frequency option
    out_mixer_c = '0';
end

time_constant = 1e-5; % [s]
demod_rate = 2e3;
% create a base configuration: disable all outputs, demods and scopes
ziDAQ('setDouble',['/' device '/demods/*/rate'], 0.0);
ziDAQ('setInt',['/' device '/demods/*/trigger'], 0);
ziDAQ('setInt',['/' device '/sigouts/*/enables/*'], 0);
if strfind(devtype,'UHF')
% if the device is a UHF additionally disable all demodulators
    ziDAQ('setInt',['/' device '/demods/*/enable'], 0);
    ziDAQ('setInt',['/' device '/scopes/*/enable'], 0);
elseif strfind(devtype,'HF2')
    ziDAQ('setInt',['/' device '/scopes/*/trigchannel'],-1)
end

% configure the device ready for this experiment
ziDAQ('setInt',['/' device '/sigins/' in_c '/imp50'], 0);
ziDAQ('setInt',['/' device '/sigins/' in_c '/ac'], 1);
ziDAQ('setDouble',['/' device '/sigins/' in_c '/range'], 1.2);
ziDAQ('setInt',['/' device '/sigouts/' out_c '/on'], 1);
ziDAQ('setDouble',['/' device '/sigouts/' out_c '/range'], 1);
ziDAQ('setDouble',['/' device '/sigouts/' out_c '/amplitudes/*'], 0);
ziDAQ('setDouble',['/' device '/sigouts/' out_c '/amplitudes/' out_mixer_c], sigout_amplitude);
ziDAQ('setDouble',['/' device '/sigouts/' out_c '/enables/' out_mixer_c], 1);
if strfind(devtype,'HF2')
    ziDAQ('setInt',['/' device '/sigins/' in_c '/diff'], 0);
    ziDAQ('setInt',['/' device '/sigouts/' out_c '/add'], 0);
end
ziDAQ('setDouble',['/' device '/demods/*/phaseshift'], 0);
ziDAQ('setInt',['/' device '/demods/*/order'], 8);
ziDAQ('setDouble',['/' device '/demods/' demod_c '/rate'], demod_rate);
ziDAQ('setInt',['/' device '/demods/' demod_c '/harmonic'], 1);
if strfind(devtype,'UHF')
    ziDAQ('setInt',['/' device '/demods/' demod_c '/enable'], 1);
end
if strfind(options,'MF')
% HF2IS and HF2LI multi-frequency option do not support the node oscselect.
    ziDAQ('setInt',['/' device '/demods/*/oscselect'], str2double(osc_c));
    ziDAQ('setInt',['/' device '/demods/*/adcselect'], str2double(in_c));
end
ziDAQ('setDouble',['/' device '/demods/*/timeconstant'], time_constant);
ziDAQ('setDouble',['/' device '/oscs/' osc_c '/freq'], freq); % [Hz]

% Unsubscribe from any streaming data
ziDAQ('unsubscribe','*');

% Pause to get a settled demodulator filter
pause(10*time_constant);

% Perform a global synchronisation between the device and the data server:
% Ensure that the settings have taken effect on the device before issuing the
% ``getSample`` command. Note: ``sync`` must be issued after waiting for the
% demodulator filter to settle above.
ziDAQ('sync');

% obtain a single demodulator sample
sample = ziDAQ('getSample',['/' device '/demods/0/sample']);
r = sqrt( sample.x^2 + sample.y^2 );
fprintf('Measured RMS amplitude is %fV.\n',r)

end

function test_example(r,amplitude)
% TEST_EXAMPLE run some simple tests on the data returned by the example
r_expected = amplitude / sqrt(2);
assert( abs(r - r_expected)/r_expected < 1, ...
        'Assertion failed: abs(r - r_expected)/r_expected < 1\n r = %.5f\n r_expected = %.5f', ...
        r, r_expected);
end

% Local variables:
% matlab-indent-level: 4
% matlab-indent-function-body: nil
% End:
