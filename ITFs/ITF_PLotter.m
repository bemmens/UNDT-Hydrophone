plot(FP15622.ITF)
hold on
plot(FP15807.ITF)
plot(FP15801.ITF)
plot(TFS564914.ITF,'--')
plot(TFS564904.ITF,'--')
hold off
ylabel('Reflectivity [%]')
xlabel('Channel')
title('ITF')
legend('FP15622','FP15807','FP15801','TFS564914','TFS564904')