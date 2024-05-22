%
% Copyright 2009-2010, Zurich Instruments AG, Switzerland
% This software is a preliminary version. Function calls and
% paramters may change without notice.
%
function [r] = multiFreqCmplx(device, n, no_frqs)
% clear ziDAQ
% ziDAQ('connect');
% device = autoDetect;
for ff = 1:1:no_frqs
    for ii = 1:1:n
        demod_c = num2str(ff-1);
        sample = ziDAQ('getSample',['/' device '/demods/' demod_c '/sample']);
        r1(ii) = sample.x + 1i.*sample.y;
    end
    frq2(ff) = sample.frequency;
    r(ff) = mean(r1); %sample.x + 1i.*sample.y;
end
clear sample
end
