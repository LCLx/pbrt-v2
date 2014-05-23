Random = [1,2,4,8,16,32,64];
LowDis = Random;
RandomMSE = [0.0326903, 0.012796, 0.00658622, 0.00344551, ...
    0.001577, 0.000929007,0.000499921];
LowDisMSE = [0.0266358, 0.00782313, 0.00182847, 0.000587319, ...
    0.000163875, 7.01577e-05, 2.10761e-05];
plot(log2(Random), log2(RandomMSE), 'r+', ...
    log2(LowDis), log2(LowDisMSE), 'g+');
xlabel('log2(samples)');
ylabel('log2(mse)');
legend('random', 'low discrepancy');
title('mse-samples: log-log scale');

