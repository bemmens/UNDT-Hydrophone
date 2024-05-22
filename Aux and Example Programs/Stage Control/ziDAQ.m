%
% Copyright 2009-2015, Zurich Instruments Ltd, Switzerland
% This software is a preliminary version. Function calls and
% parameters may change without notice.
%
% This version is linked against Matlab 7.9.0.529 (R2009b) libraries.
% In case of incompatibility with the currently used Matlab version. Send
% us the requested version (use 'ver' command).
%
% ziDAQ is a interface for communication with the ZI server.
% Usage: ziDAQ([command], [option1], [option2])
%        [command] = 'clear', 'connect', 'connectDevice',
%                    'disconnectDevice', 'finished', 'flush', 'get',
%                    'getAsEvent', 'getAuxInSample','getByte',
%                    'getDIO', 'getDouble', 'getInt',
%                    'getSample', 'listNodes', 'logOn', 'logOff',
%                    'poll', 'pollEvent', 'programRT', 'progress', 'read',
%                    'record', 'setByte', 'setDouble', 'syncSetDouble',
%                    'setInt', 'syncSetInt', 'subscribe',
%                    'sweep', 'trigger', 'unsubscribe', 'update',
%                    'zoomFFT', 'deviceSettings'
%
% Preconditions: ZI Server must be running (check task manager)
%
%            ziDAQ('connect', [host = '127.0.0.1'], [port = 8005], [apiLevel = 1]);
%                  [host] = Server host string (default is the localhost)
%                  [port] = Port number (double)
%                           Use port 8005 to connect to the HF2 Data Server
%                           Use port 8004 to connect to the UHF Data Server
%                  [apiLevel] = Compatibility mode of the API interface
%                           Use API level 1 to use code written for HF2.
%                           Higher API levels are currently only supported
%                           for UHF devices. To get full functionality for
%                           an UHF device use API level 5.
%                  To disconnect use 'clear ziDAQ'
%
%   result = ziDAQ('getConnectionAPILevel');
%                  Returns ziAPI level used for the active connection.
%
%            ziDAQ('connectDevice', [device], [interface]);
%                  [device] = Device serial string to connect (e.g. 'DEV2000')
%                  [interface] = Interface string e.g. 'USB', '1GbE', '10GbE'
%                  Connect with the data server to a specified device over the
%                  specified interface. The device must be visible to the server.
%                  If the device is already connected the call will be ignored.
%                  The function will block until the device is connected and
%                  the device is ready to use. This method is useful for UHF
%                  devices offering several communication interfaces.
%
%            ziDAQ('disconnectDevice', [device]);
%                  [device] = Device serial string of device to disconnect.
%                  This function will return immediately. The disconnection of
%                  the device may not yet finished.
%
%   result = ziDAQ('listNodes', [path], [flags(int64) = 0]);
%                  [path] = Path string
%                  [flags] = int64(0) -> ZI_LIST_NONE 0x00
%                              The default flag, returning a simple
%                              listing if the given node
%                            int64(1) -> ZI_LIST_RECURSIVE 0x01
%                              Returns the nodes recursively
%                            int64(2) -> ZI_LIST_ABSOLUTE 0x02
%                              Returns absolute paths
%                            int64(4) -> ZI_LIST_LEAFSONLY 0x04
%                              Returns only nodes that are leafs,
%                              which means the they are at the
%                              outermost level of the tree.
%                            int64(8) -> ZI_LIST_SETTINGSONLY 0x08
%                              Returns only nodes which are marked
%                              as setting
%                  Or combinations of flags might be used.
%
%   result = ziDAQ('getSample', [path]);
%                  [path] = Path string
%                  Returns a single demodulator sample (including
%                  DIO and AuxIn). For more efficient data recording
%                  use subscribe and poll functions.
%
%   result = ziDAQ('getAuxInSample', [path]);
%                  [path] = Path string
%                  Returns a single auxin sample. The auxin data
%                  is averaged in contrast to the auxin data embedded
%                  in the demodulator sample.
%
%   result = ziDAQ('getDIO', [path]);
%                  [path] = Path string
%                  Returns a single DIO sample.
%
%   result = ziDAQ('getDouble', [path]);
%                  [path] = Path string
%
%   result = ziDAQ('getInt', [path]);
%                  [path] = Path string
%
%   result = ziDAQ('getByte', [path]);
%                  [path] = Path string
%
%            ziDAQ('setDouble', [path], [value(double)]);
%                  [path] = Path string
%                  [value] = Setting value
%
%            ziDAQ('syncSetDouble', [path], [value(double)]);
%                  [path] = Path string
%                  [value] = Setting value
%
%            ziDAQ('setInt', [path], [value(int64)]);
%                  [path] = Path string
%                  [value] = Setting value
%
%            ziDAQ('syncSetInt', [path], [value(int64)]);
%                  [path] = Path string
%                  [value] = Setting value
%
%            ziDAQ('setByte', [path], [value(uint8)]);
%                  [path] = Path string
%                  [value] = Setting value
%
%            ziDAQ('subscribe', [path]);
%                  [path] = Path string
%                  Subscribe to the specified path to receive streaming data
%                  or setting data if changed over the poll or pollEvent
%                  commands.
%
%            ziDAQ('unsubscribe', [path]);
%                  [path] = Path string
%
%            ziDAQ('getAsEvent', [path]);
%                  [path] = Path string
%                  Triggers a single event on the path to return the current
%                  value. The result can be fetched with the poll or pollEvent
%                  command.
%
%            ziDAQ('update');
%                  Detect HF2 devices connected to the USB. On Windows this
%                  update is performed automatically.
%
%            ziDAQ('get', [path]);
%                  [path] = Path string
%                  Gets a structure of the node data from the specified
%                  branch. High-speed streaming nodes (e.g. /devN/demods/0/sample)
%                  are not returned. Wildcards (*) may be used, in which case
%                  read-only nodes are ignored.
%
%            ziDAQ('flush');
%                  Flush all data in the socket connection and API buffers.
%                  Call this function before a subscribe with subsequent poll
%                  to get rid of old streaming data that might still be in
%                  the buffers.
%
%            ziDAQ('echoDevice', [device]);
%                  [device] = device string e.g. 'dev100'
%                  Sends an echo command to a device and blocks until
%                  answer is received. This is useful to flush all
%                  buffers between API and device to enforce that
%                  further code is only executed after the device executed
%                  a previous command.
%
%            ziDAQ('sync');
%                  Synchronize all data path. Ensures that get and poll
%                  commands return data which was recorded after the
%                  setting changes in front of the sync command. This
%                  sync command replaces the functionality of all syncSet,
%                  flush, and echoDevice commands.
%
%            ziDAQ('programRT', [device], [filename]);
%                  [device] = device string e.g. 'dev100'
%                  [filename] = filename of RT program
%                  Writes down the RT program. To use this function
%                  the RT option must be available for the specified
%                  device.
%
%   result = ziDAQ('secondsTimeStamp', [timestamps]);
%                  [timestamps] = vector of uint64 ticks
%                  Deprecated. In order to convert timestamps to seconds divide the
%                  timestamps by the value instrument's clockbase device node,
%                   e.g., /dev99/clockbase.
%                  [Converts a timestamp vector of uint64 ticks
%                  into a double vector of timestamps in seconds (HF2 Series).]
%
% Synchronous Interface
%
%            ziDAQ('poll', [duration(double)], [timeout(int64)], [flags(uint32)]);
%                  [duration] = Recording time in [s]
%                  [timeout] = Poll timeout in [ms]
%                  [flags] = Flags that specify data polling properties
%                            Bit[0] FILL : Fill data loss holes
%                            Bit[1] ALIGN : Align data of several demodulators
%                            Bit[2] THROW : Throw if data loss is detected
%                  Records data for the specified time. This function call
%                  is blocking. Use the asynchronous interface for long
%                  recording durations.
%
%   result = ziDAQ('pollEvent', [timeout(int64)]);
%                  [timeout] = Poll timeout in [ms]
%                  Just execute a single poll command. This is a low-level
%                  function. The poll function is better suited for most
%                  cases.
%
% Asynchronous Interface
%
%   Trigger Parameters
%     trigger/buffersize      double Overwrite the buffersize [s] of the trigger
%                                    object (set when it was instantiated).
%                                    The recommended buffer size is
%                                    2*trigger/0/duration.
%     trigger/device          string The device ID to execute the software trigger,
%                                    e.g. dev123 (compulsory parameter).
%     trigger/endless         bool   Enable endless triggering 1=enable; 0=disable.
%     trigger/forcetrigger    bool   Force a trigger.
%     trigger/0/path          string The path to the demod sample to trigger on,
%                                    e.g. demods/3/sample, see also trigger/0/source
%     trigger/0/source        int    Signal that is used to trigger on.
%                                    0 = x [X_SOURCE]
%                                    1 = y [Y_SOURCE]
%                                    2 = r [R_SOURCE]
%                                    3 = angle [ANGLE_SOURCE]
%                                    4 = frequency [FREQUENCY_SOURCE]
%                                    5 = phase [PHASE_SOURCE]
%                                    6 = auxiliary input 0 [AUXIN0_SOURCE]
%                                    7 = auxiliary input 1 [AUXIN1_SOURCE]
%     trigger/0/count         int    Number of trigger edges to record.
%     trigger/0/type          int    Trigger type used. Some parameters are
%                                    only valid for special trigger types.
%                                    0 = trigger off
%                                    1 = analog edge trigger on source
%                                    2 = digital trigger mode on DIO
%                                    3 = analog pulse trigger on source
%                                    4 = analog tracking trigger on source
%     trigger/0/edge          int    Trigger edge
%                                    1 = rising edge
%                                    2 = falling edge
%                                    3 = both
%     trigger/0/findlevel     bool   Automatically find the value of trigger/0/level
%                                    based on the current signal value.
%     trigger/0/bits          int    Digital trigger condition.
%     trigger/0/bitmask       int    Bit masking for bits used for
%                                    triggering. Used for digital trigger.
%     trigger/0/delay         double Trigger frame position [s] (left side)
%                                    relative to trigger edge.
%                                    delay = 0 -> trigger edge at left border.
%                                    delay < 0 -> trigger edge inside trigger
%                                                 frame (pretrigger).
%                                    delay > 0 -> trigger edge before trigger
%                                                 frame (posttrigger).
%     trigger/0/duration      double Recording frame length [s]
%     trigger/0/level         double Trigger level voltage [V].
%     trigger/0/hysteresis    double Trigger hysteresis [V].
%     trigger/0/retrigger     int    Record more than one trigger in a
%     trigger/triggered       bool   Has the software trigger triggered? 1=Yes, 0=No
%                                    (read only).
%     trigger/0/bandwidth     double Filter bandwidth [Hz] for pulse and
%                                    tracking triggers.
%     trigger/0/holdoff/count int    Number of skipped triggers until the
%                                    next trigger is recorded again.
%     trigger/0/holdoff/time  double Hold off time [s] before the next
%                                    trigger is recorded again. A hold off
%                                    time smaller than the duration will
%                                    produce overlapped trigger frames.
%     trigger/N/hwtrigsource  int    Only available for devices that support
%                                    hardware triggering. Specify the channel
%                                    to trigger on.
%     trigger/0/pulse/min     double Minimal pulse width [s] for the pulse
%                                    trigger.
%     trigger/0/pulse/max     double Maximal pulse width [s] for the pulse
%                                    trigger.
%     trigger/filename        string This parameter is deprecated. If specified,
%                                    i.e., not empty, it enables automatic saving of
%                                    data in single trigger mode
%                                    (trigger/endless = 0).
%     trigger/savepath        string The directory where files are saved when saving
%                                    data.
%     trigger/fileformat      string The format of the file for saving data.
%                                    0=Matlab, 1=CSV.
%     trigger/historylength   bool   Maximum number of entries stored in the
%                                    measurement history.
%     triggerclearhistory     bool   Remove all records from the history list.
%
%   handle = ziDAQ('record' [duration(double)], [timeout(int64)]);
%                  [duration] = Recording time in [s]
%                  [timeout] = Poll timeout in [ms]
%                  Creates a recorder class. The thread is not yet started.
%                  Before the thread start subscribe and set command have
%                  to be called. To start the real measurement use the
%                  execute function. After that the trigger will start
%                  the recording of a frame.
%
%   result = ziDAQ('listNodes', [handle], [path], [flags(int64) = 0]);
%                  [path] = Path string
%                  [flags] = int64(0) -> ZI_LIST_NONE 0x00
%                              The default flag, returning a simple
%                              listing if the given node
%                            int64(1) -> ZI_LIST_RECURSIVE 0x01
%                              Returns the nodes recursively
%                            int64(2) -> ZI_LIST_ABSOLUTE 0x02
%                              Returns absolute paths
%                            int64(4) -> ZI_LIST_LEAFSONLY 0x04
%                              Returns only nodes that are leafs,
%                              which means the they are at the
%                              outermost level of the tree.
%                            int64(8) -> ZI_LIST_SETTINGSONLY 0x08
%                              Returns only nodes which are marked
%                              as setting
%                 Or combinations of flags might be used.
%
%            ziDAQ('subscribe', [handle], [path]);
%                  Subscribe to one or several nodes. After subscription
%                  the recording process can be started with the 'execute'
%                  command. During the recording process paths can not be
%                  subscribed or unsubscribed.
%                  [handle] = Reference to the ziDAQRecorder class.
%                  [path] = Path string of the node. Use wild card to
%                  select all. Alternatively also a list of path
%                  strings can be specified.
%
%            ziDAQ('unsubscribe', [handle], [path]);
%                  Unsubscribe from one or several nodes. During the
%                  recording process paths can not be subscribed or
%                  unsubscribed.
%                  [handle] = Reference to the ziDAQRecorder class.
%                  [path] = Path string of the node. Use wild card to
%                  select all. Alternatively also a list of path
%                  strings can be specified.
%
%            ziDAQ('get', [handle], [path]);
%                  [handle] = Reference to the ziDAQRecorder class.
%                  [path] = Path string of the node.
%
%            ziDAQ('set', [handle], [path], [value]);
%                  [handle] = Reference to the ziDAQRecorder class.
%                  [path] = Path string of the node.
%
%            ziDAQ('set', [handle], [path], [value]);
%                  [handle] = Reference to the ziDAQRecorder class.
%                  [path] = Path string of the node.
%
%            ziDAQ('set', [handle], [path], [value]);
%                  [handle] = Reference to the ziDAQRecorder class.
%                  [path] = Path string of the node.
%
%            ziDAQ('execute', [handle]);
%                  Start the recorder. After that command any trigger will
%                  start the measurement. Subscription or unsubscription
%                  is no more possible until the recording is finished.
%
%            ziDAQ('trigger', [handle]);
%                  [handle] = Handle of the recording session.
%                  Triggers the measurement recording.
%
%   result = ziDAQ('finished', [handle]);
%                  [handle] = Handle of the recording session.
%                  Returns 1 if the recording is finished, otherwise 0.
%
%   result = ziDAQ('read', [handle]);
%                  [handle] = Handle of the recording session.
%                  Transfer the recorded data to Matlab.
%
%            ziDAQ('finish', [handle]);
%                  Stop recording. The recording may be restarted by
%                  calling 'execute' again.
%
%   result = ziDAQ('progress', [handle]);
%                  Report the progress of the measurement with a number
%                  between 0 and 1.
%
%            ziDAQ('clear', [handle]);
%                  [handle] = Handle of the recording session.
%                  Stop recording.
%
% Sweep Module
%
%   Sweep Parameters
%     sweep/device           string  Device that should be used for
%                                    the parameter sweep, e.g. 'dev99'.
%     sweep/start            double  Sweep start frequency [Hz]
%     sweep/stop             double  Sweep stop frequency [Hz]
%     sweep/gridnode         string  Path of the node that should be
%                                    used for sweeping. For frequency
%                                    sweep applications this will be e.g.
%                                    'oscs/0/freq'. The device name of
%                                    the path can be omitted and is given
%                                    by sweep/device.
%     sweep/loopcount        int     Number of sweep loops (default 1)
%     sweep/endless          int     Endless sweeping (default 0)
%                                    0 = Use loopcount value
%                                    1 = Endless sweeping enabled, ignore
%                                        loopcount
%     sweep/samplecount      int     Number of samples per sweep
%     sweep/settling/time    double  Settling time before measurement is
%                                    performed, in [s]
%     sweep/settling/tc      double  Settling precision
%                                    5 ~ low precision
%                                    15 ~ medium precision
%                                    50 ~ high precision
%     sweep/settling/inaccuracy int  Demodulator filter settling inaccuracy defining
%                                    the wait time between a sweep parameter change
%                                    and recording of the next sweep point. Typical
%                                    inaccuracy values: 10m for highest sweep speed
%                                    for large signals, 100u for precise amplitude
%                                    measurements, 100n for precise noise
%                                    measurements. Depending on the order the
%                                    settling accuracy will define the number of
%                                    filter time constants the sweeper has to
%                                    wait. The maximum between this value and the
%                                    settling time is taken as wait time until the
%                                    next sweep point is recorded.
%     sweep/xmapping         int     Sweep mode
%                                    0 = linear
%                                    1 = logarithmic
%     sweep/scan             int     Scan type
%                                    0 = sequential
%                                    1 = binary
%                                    2 = bidirectional
%                                    3 = reverse
%     sweep/bandwidth        double  Fixed bandwidth [Hz]
%                                    0 = Automatic calculation (obsolete)
%     sweep/bandwidthcontrol int     Sets the bandwidth control mode (default 2)
%                                    0 = Manual (user sets bandwidth and order)
%                                    1 = Fixed (uses fixed bandwidth value)
%                                    2 = Auto (calculates best bandwidth value)
%                                        Equivalent to the obsolete bandwidth = 0
%                                        setting
%     sweep/order            int     Defines the filter roll off to use in Fixed
%                                    bandwidth selection.
%                                    Valid values are between 1 (6 dB/octave) and
%                                    8 (48 dB/octave).
%     sweep/maxbandwidth     double  Maximal bandwidth used in auto bandwidth
%                                    mode in [Hz]. The default is 1.25MHz.
%     sweep/omegasuppression double  Damping in [dB] of omega and 2omega components.
%                                    Default is 40dB in favor of sweep speed.
%                                    Use higher value for strong offset values or
%                                    3omega measurement methods.
%     sweep/averaging/tc     double  Min averaging time [tc]
%                                    0 = no averaging (see also time!)
%                                    5 ~ low precision
%                                    15 ~ medium precision
%                                    50 ~ high precision
%     sweep/averaging/sample int     Min samples to average
%                                    1 = no averaging (if averaging/tc = 0)
%     sweep/phaseunwrap      bool    Enable unwrapping of slowly changing phase
%                                    evolutions around the +/-180 degree boundary.
%     sweep/sincfilter       bool    Enables the sinc filter if the sweep frequency
%                                    is below 50 Hz. This will improve the sweep
%                                    speed at low frequencies as omega components
%                                    do not need to be suppressed by the normal
%                                    low pass filter.
%     sweep/filename          string This parameter is deprecated. If specified,
%                                    i.e. not empty, it enables automatic saving of
%                                     data in single sweep mode (sweep/endless = 0).
%     sweep/savepath          string The directory where files are located when
%                                    saving sweeper measurements.
%     sweep/fileformat        string The format of the file for saving sweeper
%                                    measurements. 0=Matlab, 1=CSV.
%     sweep/historylength     bool   Maximum number of entries stored in the
%                                    measurement history.
%     sweep/clearhistory      bool   Remove all records from the history list.
%
%     Note:
%     Settling time = max(settling.tc * tc, settling.time)
%     Averaging time = max(averaging.tc * tc, averaging.sample / sample-rate)
%
%   handle = ziDAQ('sweep', [timeout(int64)]);
%                  [timeout] = Poll timeout in [ms]
%                  Creates a sweep class. The thread is not yet started.
%                  Before the thread start subscribe and set command have
%                  to be called. To start the real measurement use the
%                  execute function.
%
%   result = ziDAQ('listNodes', [handle], [path], [flags(int64) = 0]);
%                  [path] = Path string
%                  [flags] = int64(0) -> ZI_LIST_NONE 0x00
%                              The default flag, returning a simple
%                              listing if the given node
%                            int64(1) -> ZI_LIST_RECURSIVE 0x01
%                              Returns the nodes recursively
%                            int64(2) -> ZI_LIST_ABSOLUTE 0x02
%                              Returns absolute paths
%                            int64(4) -> ZI_LIST_LEAFSONLY 0x04
%                              Returns only nodes that are leafs,
%                              which means the they are at the
%                              outermost level of the tree.
%                            int64(8) -> ZI_LIST_SETTINGSONLY 0x08
%                              Returns only nodes which are marked
%                              as setting
%                 Or combinations of flags might be used.
%
%            ziDAQ('subscribe', [handle], [path]);
%                  Subscribe to one or several nodes. After subscription
%                  the recording process can be started with the 'execute'
%                  command. During the recording process paths can not be
%                  subscribed or unsubscribed.
%                  [handle] = Reference to the ziDAQSweeper class.
%                  [path] = Path string of the node. Use wild card to
%                  select all. Alternatively also a list of path
%                  strings can be specified.
%
%            ziDAQ('unsubscribe', [handle], [path]);
%                  Unsubscribe from one or several nodes. During the
%                  recording process paths can not be subscribed or
%                  unsubscribed.
%                  [handle] = Reference to the ziDAQSweeper class.
%                  [path] = Path string of the node. Use wild card to
%                  select all. Alternatively also a list of path
%                  strings can be specified.
%
%            ziDAQ('execute', [handle]);
%                  Start the sweep. Subscription or unsubscription
%                  is no more possible until the sweep is finished.
%
%   result = ziDAQ('finished', [handle]);
%                  [handle] = Handle of the sweep session.
%                  Returns 1 if the sweep is finished, otherwise 0.
%
%   result = ziDAQ('read', [handle]);
%                  [handle] = Handle of the sweep session.
%                  Transfer the sweep data to Matlab.
%
%   result = ziDAQ('progress', [handle]);
%                  Report the progress of the measurement with a number
%                  between 0 and 1.
%
%            ziDAQ('finish', [handle]);
%                  Stop the sweep. The sweep may be restarted by
%                  calling 'execute' again.
%
%            ziDAQ('clear', [handle]);
%                  [handle] = Handle of the sweep session.
%                  Stop the current sweep.
%
%            ziDAQ('save', [handle]);
%                  Save the measured data to a file.
%                  [handle] = Handle of the sweep session.
%                  [filename] =  File in which to store the data.
%
% Zoom FFT Module
%
%   Zoom FFT Parameters
%     zoomFFT/device        string  Device that should be used for
%                                   the zoom FFT, e.g. 'dev99'.
%     zoomFFT/bit           int     Number of FFT points 2^bit
%     zoomFFT/mode          int     Zoom FFT mode
%                                   0 = Perform FFT on X+iY
%                                   1 = Perform FFT on R
%                                   2 = Perform FFT on Phase
%     zoomFFT/loopcount     int     Number of zoom FFT loops (default 1)
%     zoomFFT/endless       int     Perform endless zoom FFT (default 0)
%                                   0 = Use loopcount value
%                                   1 = Endless zoom FFT enabled, ignore
%                                       loopcount
%     zoomFFT/overlap       double  FFT overlap 0 = none, [0..1]
%     zoomFFT/settling/time double  Settling time before measurement is performed
%     zoomFFT/settling/tc   double  Settling time in time constant units before
%                                   the FFT recording is started.
%                                   5 ~ low precision
%                                   15 ~ medium precision
%                                   50 ~ high precision
%     zoomFFT/window        int     FFT window (default 1 = Hann)
%                                   0 = Rectangular
%                                   1 = Hann
%                                   2 = Hamming
%                                   3 = Blackman Harris 4 term
%     zoomFFT/absolute      bool    Shifts the frequencies so that the center
%                                   frequency becomes the demodulation frequency
%                                   rather than 0 Hz.
%
%   handle = ziDAQ('zoomFFT', [timeout(int64)]);
%                  [timeout] = Poll timeout in [ms]
%                  Creates a zoom FFT class. The thread is not yet started.
%                  Before the thread start subscribe and set command have
%                  to be called. To start the real measurement use the
%                  execute function.
%
%   result = ziDAQ('listNodes', [handle], [path], [flags(int64) = 0]);
%                  [path] = Path string
%                  [flags] = int64(0) -> ZI_LIST_NONE 0x00
%                              The default flag, returning a simple
%                              listing if the given node
%                            int64(1) -> ZI_LIST_RECURSIVE 0x01
%                              Returns the nodes recursively
%                            int64(2) -> ZI_LIST_ABSOLUTE 0x02
%                              Returns absolute paths
%                            int64(4) -> ZI_LIST_LEAFSONLY 0x04
%                              Returns only nodes that are leafs,
%                              which means the they are at the
%                              outermost level of the tree.
%                            int64(8) -> ZI_LIST_SETTINGSONLY 0x08
%                              Returns only nodes which are marked
%                              as setting
%                 Or combinations of flags might be used.
%
%            ziDAQ('subscribe', [handle], [path]);
%                  Subscribe to one or several nodes. After subscription
%                  the recording process can be started with the 'execute'
%                  command. During the recording process paths can not be
%                  subscribed or unsubscribed.
%                  [handle] = Reference to the ziDAQZoomFFT class.
%                  [path] = Path string of the node. Use wild card to
%                  select all. Alternatively also a list of path
%                  strings can be specified.
%
%            ziDAQ('unsubscribe', [handle], [path]);
%                  Unsubscribe from one or several nodes. During the
%                  recording process paths can not be subscribed or
%                  unsubscribed.
%                  [handle] = Reference to the ziDAQZoomFFT class.
%                  [path] = Path string of the node. Use wild card to
%                  select all. Alternatively also a list of path
%                  strings can be specified.
%
%            ziDAQ('execute', [handle]);
%                  Start the zoom FFT. Subscription or unsubscription
%                  is no more possible until the zoomFFT is finished.
%
%   result = ziDAQ('finished', [handle]);
%                  [handle] = Handle of the zoom FFT session.
%                  Returns 1 if the zoom FFT is finished, otherwise 0.
%
%   result = ziDAQ('read', [handle]);
%                  [handle] = Handle of the zoom FFT session.
%                  Transfer the zoomFFT data to Matlab.
%
%   result = ziDAQ('progress', [handle]);
%                  Report the progress of the measurement with a number
%                  between 0 and 1.
%
%            ziDAQ('finish', [handle]);
%                  Stop the zoomFFT. The zoom FFT may be restarted by
%                  calling 'execute' again.
%
%            ziDAQ('clear', [handle]);
%                  [handle] = Handle of the zoom FFT session.
%                  Stop the current zoom FFT.
%
% Device Settings Module
%
%   Device Settings Parameters
%     devicesettings/device         string  Device whose settings are to be
%                                           saved/loaded, e.g. 'dev99'.
%     devicesettings/path           string  Path where the settings files are to
%                                           be located. If not set, the default
%                                           settings location of the LabOne
%                                           software is used.
%     devicesettings/filename       string  The file to which the settings are to
%                                           be saved/loaded.
%     devicesettings/command        string  The save/load command to execute.
%                                           'save' = Read device settings and save
%                                                    to file.
%                                           'load' = Load settings from file and
%                                                    write to device.
%                                           'read' = Read device settings only
%                                                    (no save).
%
%   handle = ziDAQ('deviceSettings', [timeout(int64)]);
%                  [timeout] = Poll timeout in [ms]
%                  Creates a device settings class for saving/loading device
%                  settings to/from a file. Before the thread start, set the path,
%                  filename and command parameters. To run the command, use the
%                  execute function.
%
%   result = ziDAQ('listNodes', [handle], [path], [flags(int64) = 0]);
%                  [path] = Path string
%                  [flags] = int64(0) -> ZI_LIST_NONE 0x00
%                              The default flag, returning a simple
%                              listing if the given node
%                            int64(1) -> ZI_LIST_RECURSIVE 0x01
%                              Returns the nodes recursively
%                            int64(2) -> ZI_LIST_ABSOLUTE 0x02
%                              Returns absolute paths
%                            int64(4) -> ZI_LIST_LEAFSONLY 0x04
%                              Returns only nodes that are leafs,
%                              which means the they are at the
%                              outermost level of the tree.
%                            int64(8) -> ZI_LIST_SETTINGSONLY 0x08
%                              Returns only nodes which are marked
%                              as setting
%                 Or combinations of flags might be used.
%
%            ziDAQ('subscribe', [handle], [path]);
%                  Not relevant for the device settings module.
%
%            ziDAQ('unsubscribe', [handle], [path]);
%                  Not relevant for the device settings module.
%
%            ziDAQ('execute', [handle]);
%                  Execute the command.
%
%   result = ziDAQ('finished', [handle]);
%                  [handle] = Handle of the device settings session.
%                  Returns 1 if the command is finished, otherwise 0.
%
%   result = ziDAQ('read', [handle]);
%                  [handle] = Handle of the device settings session.
%                  Transfer the device settings to Matlab.
%                  Not relevant since device settings are saved to a file.
%
%   result = ziDAQ('progress', [handle]);
%                  Report the progress of the command with a number
%                  between 0 and 1.
%
%            ziDAQ('finish', [handle]);
%                  Stop the device settings module. The module may be restarted by
%                  calling 'execute' again.
%
%            ziDAQ('clear', [handle]);
%                  [handle] = Handle of the device settings session.
%                  End the current device settings thread.
%
% PLL Advisor Module
%
%   PLL Advisor Parameters
%     pllAdvisor/bode       struct  Output parameter. Contains the resulting bode
%                                   plot of the PLL simulation.
%     pllAdvisor/calculate  int     Command to calculate values. Set to 1 to start
%                                   the calculation.
%     pllAdvisor/center     double  Center frequency of the PLL oscillator. The PLL
%                                   frequency shift is relative to this center
%                                   frequency.
%     pllAdvisor/d          int     Differential gain.
%     pllAdvisor/demodbw    int     Demodulator bandwidth used for the PLL loop
%                                   filter.
%     pllAdvisor/i          double  Integral gain.
%     pllAdvisor/mode       double  Select PLL Advisor mode. Currently only one mode
%                                   (open loop) is supported.
%     pllAdvisor/order      double  Demodulator order used for the PLL loop filter.
%     pllAdvisor/p          int     Proportional gain.
%     pllAdvisor/pllbw      int     Demodulator bandwidth used for the PLL loop
%                                   filter.
%     pllAdvisor/pm         int     Output parameter. Simulated phase margin of the
%                                   PLL with the current settings. The phase margin
%                                   should be greater than 45 deg and preferably
%                                   greater than 65 deg for stable conditions.
%     pllAdvisor/pmfreq     int     Output parameter. Simulated phase margin
%                                   frequency.
%     pllAdvisor/q          int     Quality factor. Currently not used.
%     pllAdvisor/rate       int     PLL Advisor sampling rate of the PLL control
%                                   loop.
%     pllAdvisor/stable     int     Output parameter. When 1, the PLL Advisor found
%                                   a stable solution with the given settings. When
%                                   0, revise your settings and rerun the PLL
%                                   Advisor.
%     pllAdvisor/targetbw   int     Requested PLL bandwidth. Higher frequencies may
%                                   need manual tuning.
%     pllAdvisor/targetfail int     Output parameter. 1 indicates the simulated PLL
%                                   BW is smaller than the Target BW.
%
%   handle = ziDAQ('pllAdvisor', [timeout(int64)]);
%                  [timeout] = Poll timeout in [ms]
%                  Creates a PLL Advisor class for simulating the PLL in the
%                  device. Before the thread start, set the command parameters,
%                  call execute() and then set the "calculate" parameter to start
%                  the simulation.
%
%   result = ziDAQ('listNodes', [handle], [path], [flags(int64) = 0]);
%                  [path] = Path string
%                  [flags] = int64(0) -> ZI_LIST_NONE 0x00
%                              The default flag, returning a simple
%                              listing if the given node
%                            int64(1) -> ZI_LIST_RECURSIVE 0x01
%                              Returns the nodes recursively
%                            int64(2) -> ZI_LIST_ABSOLUTE 0x02
%                              Returns absolute paths
%                            int64(4) -> ZI_LIST_LEAFSONLY 0x04
%                              Returns only nodes that are leafs,
%                              which means the they are at the
%                              outermost level of the tree.
%                            int64(8) -> ZI_LIST_SETTINGSONLY 0x08
%                              Returns only nodes which are marked
%                              as setting
%                 Or combinations of flags might be used.
%
%            ziDAQ('subscribe', [handle], [path]);
%                  Subscribe to one or several nodes.
%
%            ziDAQ('unsubscribe', [handle], [path]);
%                  Unsubscribe from one or several nodes..
%
%            ziDAQ('execute', [handle]);
%                  Start the PLL Advisor.
%
%   result = ziDAQ('finished', [handle]);
%                  [handle] = Handle of the PLL Advisor session.
%                  Returns 1 if the command is finished, otherwise 0.
%
%   result = ziDAQ('read', [handle]);
%                  [handle] = Handle of the PLL Advisor session.
%                  Read pllAdvisor data. If the simulation is still ongoing only a
%                  subset of the data is returned.
%
%   result = ziDAQ('progress', [handle]);
%                  Report the progress of the command with a number
%                  between 0 and 1.
%
%            ziDAQ('finish', [handle]);
%                  Stop the PLL Advisor module.
%
%            ziDAQ('clear', [handle]);
%                  [handle] = Handle of the PLL Advisor session.
%                  End the current PLL Advisor thread.
%
% PID Advisor Module
%
%   PID Advisor Parameters
%     pidAdvisor/advancedmode      int     Disable automatic calculation of the
%                                          start and stop value.
%     pidAdvisor/auto              int     Automatic response calculation triggered
%                                          by parameter change.
%     pidAdvisor/bode              struct  Output parameter. Contains the resulting
%                                          bode plot of the PID simulation.
%     pidAdvisor/bw                double  Output parameter. Calculated system
%                                          bandwidth.
%     pidAdvisor/calculate         int     In/Out parameter. Command to calculate
%                                          values. Set to 1 to start the
%                                          calculation.
%     pidAdvisor/display/freqstart double  Start frequency for Bode plot.
%                                          For disabled advanced mode the start
%                                          value is automatically derived from the
%                                          system properties.
%     pidAdvisor/display/freqstop  double  Stop frequency for Bode plot.
%     pidAdvisor/display/timestart double  Start time for step response.
%     pidAdvisor/display/timestop  double  Stop time for step response.
%     pidAdvisor/dut/bw            double  Bandwith of the DUT (device under test).
%     pidAdvisor/dut/damping       double  Damping of the second order
%                                          low pass filter.
%     pidAdvisor/dut/delay         double  IO Delay of the feedback system
%                                          describing the earliest response for
%                                          a step change.
%     pidAdvisor/dut/fcenter       double  Resonant frequency of the of the modelled
%                                          resonator.
%     pidAdvisor/dut/gain          double  Gain of the DUT transfer function.
%     pidAdvisor/dut/q             double  quality factor of the modelled resonator.
%     pidAdvisor/dut/source        int     Type of model used for the external
%                                          device to be controlled by the PID.
%                                          source = 1: Lowpass first order
%                                          source = 2: Lowpass second order
%                                          source = 3: Resonator frequency
%                                          source = 4: Internal PLL
%                                          source = 5: VCO
%                                          source = 6: Resonator amplitude
%     pidAdvisor/impulse           struct  Output parameter. Impulse response
%                                          (not yet supported).
%     pidAdvisor/index             int     PID index for parameter detection.
%     pidAdvisor/pid/autobw        int     Adjusts the demodulator bandwidth to fit
%                                          best to the specified target bandwidth
%                                          of the full system.
%     pidAdvisor/pid/d             double  In/Out parameter. Differential gain.
%     pidAdvisor/pid/dlimittimeconstant
%                                  double  In/Out parameter. Differential filter
%                                          timeconstant.
%     pidAdvisor/pid/i             double  In/Out parameter. Integral gain.
%     pidAdvisor/pid/mode          double  Select PID Advisor mode. Mode value is
%                                          bit coded, bit 0: P, bit 1: I, bit 2: D,
%                                          bit 3: D filter limit.
%     pidAdvisor/pid/p             double  In/Out parameter. Proportional gain.
%     pidAdvisor/pid/rate          double  In/Out parameter. PID Advisor sampling
%                                          rate of the PID control loop.
%     pidAdvisor/pid/targetbw      double  PID system target bandwidth.
%     pidAdvisor/pm                double  Output parameter. Simulated phase margin
%                                          of the PID with the current settings.
%                                          The phase margin should be greater than
%                                          45 deg and preferably greater than 65 deg
%                                          for stable conditions.
%     pidAdvisor/pmfreq            double  Output parameter. Simulated phase margin
%                                          frequency.
%     pidAdvisor/stable            int     Output parameter. When 1, the PID Advisor
%                                          found a stable solution with the given
%                                          settings. When 0, revise your settings
%                                          and rerun the PID Advisor.
%     pidAdvisor/step              struct  Output parameter. Contains the resulting
%                                          step response plot of the PID simulation.
%     pidAdvisor/targetbw          double  Requested PID bandwidth. Higher
%                                          frequencies may need manual tuning.
%     pidAdvisor/targetfail        int     Output parameter. 1 indicates the
%                                          simulated PID BW is smaller than the
%                                          Target BW.
%     pidAdvisor/tf/closedloop     int     Switch the response calculation mode
%                                          between closed or open loop.
%     pidAdvisor/tf/input          int     Start point for the plant response
%                                          simulation for open or closed loops.
%     pidAdvisor/tf/output         int     End point for the plant response
%                                          simulation for open or closed loops.
%     pidAdvisor/tune              int     Optimize the PID parameters so that
%                                          the noise of the closed-loop
%                                          system gets minimized.
%
%   handle = ziDAQ('pidAdvisor', [timeout(int64)]);
%                  [timeout] = Poll timeout in [ms]
%                  Creates a PID Advisor class for simulating the PID in the
%                  device. Before the thread start, set the command parameters,
%                  call execute() and then set the "calculate" parameter to start
%                  the simulation.
%
%   result = ziDAQ('listNodes', [handle], [path], [flags(int64) = 0]);
%                  [handle] = Handle of the PID Advisor session.
%                  [path] = Path string
%                  [flags] = int64(0) -> ZI_LIST_NONE 0x00
%                              The default flag, returning a simple
%                              listing if the given node
%                            int64(1) -> ZI_LIST_RECURSIVE 0x01
%                              Returns the nodes recursively
%                            int64(2) -> ZI_LIST_ABSOLUTE 0x02
%                              Returns absolute paths
%                            int64(4) -> ZI_LIST_LEAFSONLY 0x04
%                              Returns only nodes that are leafs,
%                              which means the they are at the
%                              outermost level of the tree.
%                            int64(8) -> ZI_LIST_SETTINGSONLY 0x08
%                              Returns only nodes which are marked
%                              as settings
%                 Or combinations of flags might be used.
%
%            ziDAQ('subscribe', [handle], [path]);
%                  [handle] = Handle of the PID Advisor session.
%                  Subscribe to one or several nodes.
%
%            ziDAQ('unsubscribe', [handle], [path]);
%                  [handle] = Handle of the PID Advisor session.
%                  Unsubscribe from one or several nodes..
%
%            ziDAQ('get', [handle], [path]);
%                  [handle] = Handle of the PID Advisor session.
%                  [path] = Path string of the node.
%
%            ziDAQ('execute', [handle]);
%                  [handle] = Handle of the PID Advisor session.
%                  Starts the pidAdvisor if not yet running.
%
%            ziDAQ('trigger', [handle]);
%                  Not applicable to this module.
%
%   result = ziDAQ('finished', [handle]);
%                  [handle] = Handle of the PID Advisor session.
%                  Returns 1 if the command is finished, otherwise 0.
%
%   result = ziDAQ('read', [handle]);
%                  [handle] = Handle of the PID Advisor session.
%                  Read pidAdvisor data. If the simulation is still ongoing only a
%                  subset of the data is returned.
%
%   result = ziDAQ('progress', [handle]);
%                  [handle] = Handle of the PID Advisor session.
%                  Report the progress of the command with a number
%                  between 0 and 1.
%
%            ziDAQ('finish', [handle]);
%                  [handle] = Handle of the PID Advisor session.
%                  Stop the PID Advisor module.
%
%            ziDAQ('clear', [handle]);
%                  [handle] = Handle of the PID Advisor session.
%                  End the current PID Advisor thread.
%
%            ziDAQ('save', [handle]);
%                  Save the measured data to a file.
%                  [handle] = Handle of the PID Advisor session.
%                  [filename] =  File name string (without extension)..
%
% Debugging Functions
%
%            ziDAQ('setDebugLevel', [debuglevel]);
%                  [debuglevel] = Debug level (trace:0, info:1, debug:2, warning:3,
%                  error:4, fatal:5, status:6).
%                  Enables debug log and sets the debug level.
%
%            ziDAQ('writeDebugLog', [severity], [message]);
%                  [severity] = Severity (trace:0, info:1, debug:2, warning:3,
%                  error:4, fatal:5, status:6).
%                  [message] = Message to output to the log.
%                  Outputs message to the debug log (if enabled).
%
%            ziDAQ('logOn', [flags], [filename], [style]);
%                  Flags = LOG_NONE:             0x00000000
%                          LOG_SET_DOUBLE:       0x00000001
%                          LOG_SET_INT:          0x00000002
%                          LOG_SET_BYTE:         0x00000004
%                          LOG_SYNC_SET_DOUBLE:  0x00000010
%                          LOG_SYNC_SET_INT:     0x00000020
%                          LOG_SYNC_SET_BYTE:    0x00000040
%                          LOG_GET_DOUBLE:       0x00000100
%                          LOG_GET_INT:          0x00000200
%                          LOG_GET_BYTE:         0x00000400
%                          LOG_GET_DEMOD:        0x00001000
%                          LOG_GET_DIO:          0x00002000
%                          LOG_GET_AUXIN:        0x00004000
%                          LOG_LISTNODES:        0x00010000
%                          LOG_SUBSCRIBE:        0x00020000
%                          LOG_UNSUBSCRIBE:      0x00040000
%                          LOG_GET_AS_EVENT:     0x00080000
%                          LOG_UPDATE:           0x00100000
%                          LOG_POLL_EVENT:       0x00200000
%                          LOG_POLL:             0x00400000
%                          LOG_ALL :             0xffffffff
%                  [filename] = Log file name
%                  [style] = LOG_STYLE_TELNET: 0 (default)
%                            LOG_STYLE_MATLAB: 1
%                            LOG_STYLE_PYTHON: 2
%                  Log all messages sent to the ziServer. This is useful
%                  for debugging.
%
%            ziDAQ('logOff');
%                  Turn of message logging.
%
