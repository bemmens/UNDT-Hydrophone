%
% Copyright 2009-2010, Zurich Instruments AG, Switzerland
% This software is a preliminary version. Function calls and
% paramters may change without notice.
%
function device = autoDetect
nodes = lower(ziDAQ('listNodes', '/'));
dutIndex = strmatch('dev', nodes);
if length(dutIndex) > 1
  error('autoDetect does only support a single device configuration.');
elseif isempty(dutIndex)
  error('No DUT found. Make sure that the USB cable is connected to the host and the device is turned on.');
end
% Found only one device -> selection valid.
device = lower(nodes{dutIndex});
fprintf('Will perform measurement for device %s ...\n', device)