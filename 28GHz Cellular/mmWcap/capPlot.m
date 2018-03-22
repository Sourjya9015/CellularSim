tests = {'txPowDLSINR', ... % 1
    'txPowSE', ...     % 2
    'bfCDF',...        % 3
    'sinrTxUL',...     % 4
    'sinrTxDL',...     % 5
    'rateTxUL',...     % 6
    'rateTxDL',...     % 7
    'rateBWDL',...     % 8
    'rateBWUL',...     % 9
    'rateMAUL',...     % 10
    'seBW',...         % 11
    'effSINRTx',...    % 12
    'inrCDF', ...      % 13
    'rateUserUL', ...  % 14
    'rateUserDL', ...   % 15
    'rateStream', ... % 16
    'sinrFreqUL', ...   % 17
    'sinrFreqDL', ...   % 18
    'rateFreqUL', ...   % 19
    'rateFreqDL', ...   % 20
    'rateBlkUL', ...   % 21
    'rateBlkDL', ...   % 22
    };
plotNames = tests([21,22]);        % select some above
useDiskData = true;         % make false to use variables in the workspace
fileSuffix = '';

figwidth = 7;
figheight = 5;

for plotind = 1 : length(plotNames);
    plotName = plotNames{plotind};
    if useDiskData
        switch plotName
            case {'txPowDLSINR','txPowSE','sinrTxUL','sinrTxDL','rateTxUL',...
                    'rateTxDL','effSINRTx','inrCDF'}
                cmd = 'load data/TxPow';
            case {'rateBW','seBW'}
                cmd = 'load data/BW';
            case 'rateMAUL'
                cmd = 'load data/powerMA';
            case {'rateUserUL','rateUserDL'}
                cmd = 'load data/user';
            case 'rateStream'
                cmd = 'load data/nstream';
            case {'sinrFreqUL','sinrFreqDL', 'rateFreqUL', 'rateFreqDL'}
                cmd = 'load data/freq';
            case {'rateBlkUL', 'rateBlkDL'}
                cmd = 'load data/block';
                
        end
        cmd = sprintf('%s%s', cmd, fileSuffix);
        eval(cmd);
    end
    
    f = figure;
    switch plotName
        
        case 'txPowDLSINR'
            % load data/TxPow;
            ntest = length(sinrTestDL);
            legSNR = cell(ntest,1);
            colors = (0:1/(ntest-1):1)';
            
            for itest = 1 : ntest
                nue = length(sinrTestDL{itest});
                p = (1:nue)/nue;
                plot(sinrTestDL{itest}, p, '-','LineWidth',2, ...
                    'Color', [colors(itest) 1-colors(itest) 1]);
                hold on;
                legSNR{itest} = sprintf('%d dBm', paramTest{itest}{2});
            end
            set(gca,'FontSize',16);
            xlabel('DL SINR (dB)');
            ylabel('Cummulative prob');
            grid on;
            axis([-30 25 0 1]);
            hold off;
            legend(legSNR);
            
        case 'txPowSE'
            % load data/TxPow;
            dir = 'DL';
            if strcmp(dir, 'UL')
                ind = 4;
                seTest = seTestUL;
                rateTest = rateTestUL;
            else
                ind = 2;
                seTest = seTestDL;
                rateTest = rateTestDL;
            end
            
            ntest = length(seTest);
            leg = cell(ntest,1);
            txPow = zeros(ntest,1);
            for itest = 1:ntest
                txPow(itest) = paramTest{itest}{ind};
                leg{itest} = sprintf('%d dBm', txPow(itest));
            end
            colors = (0:1/(ntest-1):1)';
            
            subplot(1,2,1);
            
            for itest = 1 : ntest
                nue = length(rateTest{itest});
                p = (1:nue)/nue;
                semilogx(rateTest{itest}, p, '-','LineWidth',2, ...
                    'Color', [colors(itest) 1-colors(itest) 1]);
                hold on;
            end
            set(gca,'FontSize',16);
            xlabel('Rate (Mbps)');
            ylabel('Cummulative prob');
            grid on;
            axis([0.1 1e3 0 1]);
            hold off;
            legend(leg);
            
            subplot(1,2,2);
            plot(txPow, seTest, '-o','LineWidth',2);
            set(gca,'FontSize',16);
            xlabel('TX pow (dBm)');
            ylabel('Spectral efficiency (bps/Hz/cell)');
            grid on;
            axis([min(txPow) max(txPow) 0 3]);
            
            if (strcmp(dir,'UL'))
                print -dpng plots/txPowUL;
            else
                print -dpng plots/txPowDL;
            end
            
        case 'bfCDF'
            load data/bfGain4X4.mat
            hold on;
            grid on;
            box on;
            [len, ~]=size(gainMat);
            gainMat = 10*log10(gainMat);
            plot(gainMat(:,1),(1:len)/len,'LineWidth',2,'color',[0 0 1]);
            plot(gainMat(:,2),(1:len)/len,'--','LineWidth',2,'color',[0 0 1]);
            legStr = {'Ser 4x4', 'Int 4x4'};
            load data/bfGain4X8.mat
            gainMat = 10*log10(gainMat);
            plot(gainMat(:,1),(1:len)/len,'LineWidth',2,'color',[0 1 0]);
            plot(gainMat(:,2),(1:len)/len,'--','LineWidth',2,'color',[0 1 0]);
            legStr = [legStr, {'Ser 4x8','Int 4x8'}];
            load data/bfGain8X4.mat
            gainMat = 10*log10(gainMat);
            plot(gainMat(:,1),(1:len)/len,'LineWidth',2,'color',[1 0 0]);
            plot(gainMat(:,2),(1:len)/len,'--','LineWidth',2,'color',[1 0 0]);
            legStr = [legStr, {'Ser 8x4','Int 8x4'}];
            load data/bfGain8X8.mat
            gainMat = 10*log10(gainMat);
            plot(gainMat(:,1),(1:len)/len,'LineWidth',2,'color',[1 0 1]);
            plot(gainMat(:,2),(1:len)/len,'--','LineWidth',2,'color',[1 0 1]);
            legStr = [legStr, {'Ser 8x8','Int 8x8'}];
            set(gca,'FontSize',16);
            xlabel('Beamforming Gain (dB)');
            ylabel('Empirical probability');
            title('Beamforming gain CDF');
            legend(legStr,'Location','NorthWest');
            axis([-15 20 0 1]);
            
        case 'sinrTxUL'
            % load data/TxPow.mat
            ntest = length(sinrTestUL);
            colors = (0:1/(ntest-1):1)';
            legStr = cell(ntest,1);
            power = 10:5:30;
            styles = {'-','--','-',':','-.'};
            for i = 1:ntest
                linewid = 2;   if i==3, linewid=3; end
                len = length(sinrTestUL{i});
                plot(sinrTestUL{i}, (1:len)/len, styles{i}, 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                legStr{i} = sprintf('%d dBm', power(i));
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlim([-20 30])
            xlabel('SINR (dB)');
            ylabel('Empirical probability');
            title('Uplink SINR CDF');
            legend(legStr,'Location','best');
            
            
        case 'sinrTxDL'
            % load data/TxPow.mat
            ntest = length(sinrTestDL);
            colors = (0:1/(ntest-1):1)';
            legStr = cell(ntest,1);
            power = 20:5:40;
            styles = {'-','--','-',':','-.'};
            for i = 1:ntest
                linewid = 2;   if i==3, linewid=3; end
                len = length(sinrTestDL{i});
                plot(sinrTestDL{i}, (1:len)/len, styles{i}, 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                legStr{i} = sprintf('%d dBm', power(i));
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlim([-20 30])
            xlabel('SINR (dB)');
            ylabel('Empirical probability');
            title('Downlink SINR CDF');
            legend(legStr,'Location','best');
            
            
        case 'rateTxUL'
            % load data/TxPow.mat
            ntest = length(rateTestUL);
            colors = (0:1/(ntest-1):1)';
            legStr = cell(ntest,1);
            power = 10:5:30;
            styles = {'-','--','-',':','-.'};
            for i = 1:ntest
                linewid = 2;   if i==3, linewid=3; end
                len = length(rateTestUL{i});
                semilogx(rateTestUL{i}, (1:len)/len, styles{i}, 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                legStr{i} = sprintf('%d dBm', power(i));
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            axis([1 1e3 0 1]);
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Uplink rate CDF');
            legend(legStr,'Location','best');
            
            
        case 'rateTxDL'
            % load data/TxPow.mat
            ntest = length(rateTestDL);
            colors = (0:1/(ntest-1):1)';
            legStr = cell(ntest,1);
            power = 20:5:40;
            styles = {'-','--','-',':','-.'};
            for i = 1:ntest
                linewid = 2;   if i==3, linewid=3; end
                len = length(rateTestDL{i});
                semilogx(rateTestDL{i}, (1:len)/len, styles{i}, 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                legStr{i} = sprintf('%d dBm', power(i));
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlim([1, 1e3])
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Downlink rate CDF');
            legend(legStr,'Location','best');
            
            
        case 'rateBWDL'
            % load data/BW.mat
            ntest = length(rateTestDL);
            colors = (0:1/(ntest-1):1)';
            styles = {'--','-.',':','-'};
            for i = 1:ntest
                len = length(rateTestDL{i});
                semilogx(rateTestDL{i}, (1:len)/len, styles{i}, 'LineWidth',2,'color', ...
                    [colors(i) 1-colors(i) 1]);
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlim([1, 1e3])
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Downlink rate CDF');
            legend(legStr,'Location','Northwest');
        case 'rateBWUL'
            ntest = length(rateTestUL);
            colors = (0:1/(ntest-1):1)';
            for i = 1:ntest
                len = length(rateTestUL{i});
                semilogx(rateTestUL{i}, (1:len)/len, styles{i},'LineWidth',2,'color', ...
                    [colors(i) 1-colors(i) 1]);
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlim([1, 1e3])
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Uplink rate CDF');
            legend(legStr,'Location','Northwest');
            
            
        case 'rateMAUL'
            % load data/powerMA.mat
            ntest = length(rateTestUL);
            colors = (0:1/(ntest-1):1)';
            styles = {'--','-.',':','-'};
            for i = 1:ntest
                len = length(rateTestUL{i});
                semilogx(rateTestUL{i}, (1:len)/len, styles{i},'LineWidth',2,'color', ...
                    [colors(i) 1-colors(i) 1]);
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            axis([0.1 1e3 0 1]);
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Uplink rate CDF');
            legend(legStr,'Location','Northwest');
            
            
        case 'seBW'
            % load data/BW.mat
            hold on;
            xval = [];
            for i = 1:length(paramTest)
                xval = [xval, paramTest{i}{1,2}];
            end
            plot(xval, seTestDL,'LineWidth',2,'color', [0 0 1]);
            plot(xval, seTestUL,'--','LineWidth',2,'color', [0 1 0]);
            grid on;
            box on;
            set(gca,'FontSize',16);
            axis([100,1000,1.8,3]);
            xlabel('Bandwidth (MHz)');
            ylabel('Spectral efficiency (bps/Hz/cell)');
            legend('Downlink','Uplink','Location','best');
            
        case 'effSINRTx'
            % load data/TxPow.mat
            hold on;
            xval = [];
            for i = 1:length(paramTest)
                xval = [xval, paramTest{i}{1,2}];
            end
            plot(xval, 2.^(seTestDL)-1,'LineWidth',2,'color', [0 0 1]);
            plot(xval-10, 2.^(seTestUL)-1,'--','LineWidth',2,'color', [0 1 0]);
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlabel('TX Power (dBm)');
            ylabel('Cell effective SINR');
            legend('Downlink','Uplink','Location','best');
            
            
        case 'inrCDF'
            %load data/TxPow.mat
            hold on;
            len = length(inrTestDL{3});
            plot(inrTestDL{3}, (1:len)/len,'LineWidth',2,'color', [0 0 1]);
            len = length(inrTestUL{3});
            plot(inrTestUL{3}, (1:len)/len,'--','LineWidth',2,'color', [0 1 0]);
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlim([-15, 5]);
            xlabel('INR (dB)');
            ylabel('Empirical probability');
            title('INR CDF')
            legend('Downlink 30 dBm','Uplink 20 dBm','Location','best');
            
        case 'rateUserUL'
            ntest = length(rateTestUL);
            colors = (0:1/(ntest/2-1):1)';
            legStr = cell(ntest,1);
            for i = 1:ntest/2
                linewid = 2;
                len = length(rateTestUL{i});
                semilogx(rateTestUL{i}, (1:len)/len, '--', 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                legStr{2*i-1} = sprintf('%d UEs/pico, %s', paramTest{i}{1,2},paramTest{i}{1,4});
                hold on;
                semilogx(rateTestUL{i+ntest/2}, (1:len)/len, '-', 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                legStr{2*i} = sprintf('%d UEs/pico, %s', paramTest{i+ntest/2}{1,2},paramTest{i+ntest/2}{1,4});
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            axis([1 1e3 0 1]);
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Uplink rate CDF');
            legend(legStr,'Location','best');
            
        case 'rateUserDL'
            ntest = length(rateTestDL);
            colors = (0:1/(ntest/2-1):1)';
            legStr = cell(ntest,1);
            for i = 1:ntest/2
                linewid = 2;
                len = length(rateTestDL{i});
                semilogx(rateTestDL{i}, (1:len)/len, '--', 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                legStr{2*i-1} = sprintf('%d UEs/pico, %s', paramTest{i}{1,2},paramTest{i}{1,4});
                hold on;
                semilogx(rateTestDL{i+ntest/2}, (1:len)/len, '-', 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                legStr{2*i} = sprintf('%d UEs/pico, %s', paramTest{i+ntest/2}{1,2},paramTest{i+ntest/2}{1,4});
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            axis([1 1e3 0 1]);
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Downlink rate CDF');
            legend(legStr,'Location','best');
            
        case 'rateStream'
            % load data/TxPow.mat
            dir = 'DL';
            if strcmp(dir, 'UL')
                titleStr = 'Uplink rate';
                rateTest = rateTestUL;
                plotName = 'rateStreamUL';
            else
                titleStr = 'Downlink rate';
                rateTest = rateTestDL;
                plotName = 'rateStreamDL';
            end
            ntest = length(rateTest);
            colors = (0:1/(ntest-1):1)';
            legStr = cell(ntest,1);
            for i = [1:ntest]
                nstr = paramTest{i}{2};
                linewid = 2;
                len = length(rateTest{i});
                semilogx(rateTest{i}, (1:len)/len, '-', 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                if (nstr == 0)
                    legStr{i} = sprintf('nstream = inf');
                else
                    legStr{i} = sprintf('nstream = %d', nstr);
                end
                hold on;
            end
            grid on;
            hold off;
            box on;
            set(gca,'FontSize',16);
            axis([1 1e3 0 1]);
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title(titleStr);
            legend(legStr,'Location','best');
            
        case 'sinrFreqUL'
            % load data/TxPow.mat
            ntest = length(sinrTestUL);
            colors = (0:1/(ntest-1):1)';
            legStr = cell(ntest,1);
            styles = {'-','--',':','-.','-'};
            for i = 1:ntest
                linewid = 2;
                len = length(sinrTestUL{i});
                plot(sinrTestUL{i}, (1:len)/len, styles{i}, 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                fGHz = paramTest{i}{2};
                nant1 = round(sqrt(paramTest{i}{4}));
                
                legStr{i} = sprintf('%d GHz, UE %dx%d', fGHz, nant1, nant1);
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlim([-10 30])
            xlabel('SINR (dB)');
            ylabel('Empirical probability');
            title('Uplink SINR CDF');
            legend(legStr,'Location','NorthWest');
            
            
        case 'sinrFreqDL'
            % load data/TxPow.mat
            ntest = length(sinrTestDL);
            colors = (0:1/(ntest-1):1)';
            legStr = cell(ntest,1);
            styles = {'-','--',':','-.','-'};
            for i = 1:ntest
                linewid = 2;
                len = length(sinrTestDL{i});
                plot(sinrTestDL{i}, (1:len)/len, styles{i}, 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                fGHz = paramTest{i}{2};
                nant1 = round(sqrt(paramTest{i}{4}));
                
                legStr{i} = sprintf('%d GHz, UE %dx%d', fGHz, nant1, nant1);
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlim([-10 30])
            xlabel('SINR (dB)');
            ylabel('Empirical probability');
            title('Downlink SINR CDF');
            legend(legStr,'Location','NorthWest');
            
        case 'rateFreqUL'
            % load data/TxPow.mat
            ntest = length(rateTestUL);
            colors = (0:1/(ntest-1):1)';
            legStr = cell(ntest,1);
            styles = {'-','--',':','-.','-'};
            for i = 1:ntest
                linewid = 2;
                len = length(sinrTestDL{i});
                semilogx(rateTestUL{i}, (1:len)/len, styles{i}, 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                fGHz = paramTest{i}{2};
                nant1 = round(sqrt(paramTest{i}{4}));
                
                legStr{i} = sprintf('%d GHz, UE %dx%d', fGHz, nant1, nant1);
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            axis([10 1e3 0 1]);
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Uplink rate CDF');
            legend(legStr,'Location','NorthWest');
            
            
        case 'rateFreqDL'
            % load data/TxPow.mat
            ntest = length(rateTestDL);
            colors = (0:1/(ntest-1):1)';
            legStr = cell(ntest,1);
            styles = {'-','--',':','-.','-'};
            for i = 1:ntest
                linewid = 2;
                len = length(sinrTestDL{i});
                semilogx(rateTestDL{i}, (1:len)/len, styles{i}, 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                fGHz = paramTest{i}{2};
                nant1 = round(sqrt(paramTest{i}{4}));
                
                legStr{i} = sprintf('%d GHz, UE %dx%d', fGHz, nant1, nant1);
                hold on;
            end
            grid on;
            box on;
            set(gca,'FontSize',16);
            axis([10 1e3 0 1]);
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Downlink rate CDF');
            legend(legStr,'Location','NorthWest');
            
        case 'rateBlkUL'
            % load data/TxPow.mat
            ntest = length(rateTestDL);
            colors = (0:1/(ntest-1):1)';
            %legStr = {'Threshold blocking, T=175m', 'Threshold blocking, T=300m', 'Soft blocking'};
            legStr = cell(1,ntest);
            styles = {'-','--',':','-.','-'};
            for i = 1:ntest
                linewid = 2;
                len = length(sinrTestDL{i});
                semilogx(rateTestUL{i}, (1:len)/len, styles{i}, 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                hold on;
                switch paramTest{i}{2}
                    case 'curvefit'
                        linkstate = 'Hybrid';
                    case 'curvefitNLOS'
                        linkstate = 'NLOS + Outage';
                    otherwise
                        linkstate = 'Other';
                end     
                legStr{i} = sprintf('%s, d_{shift} = %d m',linkstate,-paramTest{i}{4});
            end
            hold off;
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlim([10, 1e3])
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Uplink rate CDF');
            legend(legStr,'Location','NorthWest');
            
        case 'rateBlkDL'
            % load data/TxPow.mat
            ntest = length(rateTestDL);
            colors = (0:1/(ntest-1):1)';
            %legStr = {'Threshold blocking, T=175m', 'Threshold blocking, T=300m', 'Soft blocking'};
            legStr = cell(1,ntest);
            styles = {'-','--',':','-.','-'};
            for i = 1:ntest
                linewid = 2;
                len = length(sinrTestDL{i});
                semilogx(rateTestDL{i}, (1:len)/len, styles{i}, 'LineWidth', linewid,'color', ...
                    [colors(i) 1-colors(i) 1]);
                hold on;
                switch paramTest{i}{2}
                    case 'curvefit'
                        linkstate = 'Hybrid';
                    case 'curvefitNLOS'
                        linkstate = 'NLOS + Outage';
                    otherwise
                        linkstate = 'Other';
                end     
                legStr{i} = sprintf('%s, d_{shift} = %d m',linkstate,-paramTest{i}{4});
            end
            hold off;
            grid on;
            box on;
            set(gca,'FontSize',16);
            xlim([10, 1e3])
            xlabel('Rate (Mbps)');
            ylabel('Empirical probability');
            title('Downlink rate CDF');
            legend(legStr,'Location','NorthWest');
            
    end
    
    set(f,  'Units','inches',...
        'Position',[2 2 figwidth figheight],...
        'PaperSize',[figwidth figheight],...
        'PaperPositionMode','auto',...
        'InvertHardcopy', 'on',...
        'Renderer','painters'...
        );
    
    saveas(f,['plots/' plotName],'pdf');
    saveas(f,['plots/' plotName],'png');
    saveas(f,['plots/' plotName],'eps');
    close(f)
    
    %print (f, ['Figures/' plotName '.pdf'], '-dpdf');
    %print (f, ['Figures/', plotName, '.png'], '-dpng');
end