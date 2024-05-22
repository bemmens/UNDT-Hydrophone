function [TestArray,frq2] = MahinSingleMeasurement(TestArray,TestNo)

amp = 200e-3;
freq = 5e6;
ave = 5;

clear ziDAQ
% Connect to the Lock-in Amplifier

ziDAQ('connect');
device = autoDetect;

frq2 = example_connect_config_rob(device, freq, amp, 1);

[r] = multiFreqCmplx(device, ave, 1);

    if   TestNo == 1;
         TestArray = r;
    else TestArray = horzcat(TestArray,r);

    end

end


