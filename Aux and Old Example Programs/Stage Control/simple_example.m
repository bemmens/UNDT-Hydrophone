%
% Copyright 2009-2010, Zurich Instruments AG, Switzerland
% This software is a preliminary version. Function calls and
% paramters may change without notice.
%
function simple_example
clear ziDAQ
ziDAQ('connect');
device = autoDetect;
sample = ziDAQ('getSample', ['/' device '/demods/0/sample']);
r = sqrt(sample.x.^2 + sample.y.^2);
fprintf('Measured rms amplitude %gV\n', r);