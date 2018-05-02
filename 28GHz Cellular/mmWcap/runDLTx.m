clc;

nTests = 5;
servData = cell(nTests,1);
latData = cell(nTests,1);

multacsDL = 'fdma';
bsNumStreams = 1;
capSimDLDataTx;
servData{1} = probLat;
latData{1} = latency;

multacsDL = 'tdma';
bsNumStreams = 1;
capSimDLDataTx;
servData{2} = probLat;
latData{2} = latency;

multacsDL = 'sdma';
bsNumStreams = 2;
capSimDLDataTx;
servData{3} = probLat;
latData{3} = latency;

multacsDL = 'sdma';
bsNumStreams = 4;
capSimDLDataTx;
servData{4} = probLat;
latData{4} = latency;

multacsDL = 'sdma';
bsNumStreams = 16;
capSimDLDataTx;
servData{5} = probLat;
latData{5} = latency;

%figure(1); legend('OFDMA','TDMA, K=1','SDMA, K=2', 'SDMA, K=4', 'SDMA, K=16');
%figure(3); legend('OFDMA','TDMA, K=1','SDMA, K=2', 'SDMA, K=4', 'SDMA, K=16');

figure(4);
%OFDMA
semilogx(latData{1},servData{1}, '-', 'Linewidth',2);
hold on;

for iCnt = 2:nTests
    semilogx(latData{iCnt},servData{iCnt}, '--', 'Linewidth',2);
end
grid on;
axis([1e-4 1e2 0.5 1]);
set(gca,'FontSize',16);
xlabel('Avg. waiting time (ms)');
ylabel('Cummulative prob');
legend('OFDMA','TDMA, K=1','SDMA, K=2', 'SDMA, K=4', 'SDMA, K=16');
        