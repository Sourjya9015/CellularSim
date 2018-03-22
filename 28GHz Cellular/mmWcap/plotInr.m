f = figure;
figwidth = 7;
figheight = 5;

DLINR = inrDLs(:);
len = length(DLINR);
plot(sort(DLINR), (1:len)/len,'LineWidth',2,'color', [0 0 1]); hold on;
ULINR = inrULs(:);
len = length(ULINR);
plot(sort(ULINR), (1:len)/len,'--','LineWidth',2,'color', [0 1 0]);
grid on;
box on;
set(gca,'FontSize',16);
xlabel('INR (dB)');
ylabel('Empirical probability');
title('INR CDF')
legend('Downlink 30 dBm','Uplink 20 dBm','Location','NorthWest');
xlim([-50, 50]);
plotName = 'inrCdf';

set(f,  'Units','inches',...
    'Position',[2 2 figwidth figheight],...
    'PaperSize',[figwidth figheight],...
    'PaperPositionMode','auto',...
    'InvertHardcopy', 'on',...
    'Renderer','painters'...
    );

saveas(f,['plots/' plotName],'pdf');
saveas(f,['plots/' plotName],'png');
close(f)