% capTest:  DL and UL batch test script
% -------------------------------------

% Parameters to test
% The parameters are listed as a cell of cells.
% paramTest{itest} is the cell of parameter-value pairs for test itest

calcUL = true;      % calculate UL capacity
calcDL = true;      % calculate DL capacity
plotAvgSE = false;
savedat = true;

tests = {'antenna' ...      % 1
    'powerMA' ...           % 2
    'sectorization' ...     % 3
    'user' ...              % 4
    'TxPow' ...             % 5
    'BW' ...                % 6
    'ISD', ...              % 7
    'sect', ...             % 8
    'LOS',...               % 9
    'nstream', ...          % 10
    'freq', ...             % 11
    'block', ...            % 12
    };
testName = tests{4};        % select one above
fileSuffix = '';

switch testName
    
    case 'antenna'
        % Antenna test
        paramTest = { ...
            {'nantBS', 16, 'nantUE', 16}, ...
            {'nantBS', 64, 'nantUE', 16}, ...
            {'nantBS', 16, 'nantUE', 64}, ...
            {'nantBS', 64, 'nantUE', 64} };
        legFmtStr = 'nant=%d x %d'; % Format string for legend
        
    case 'powerMA'
        % Power test
        paramTest = { ...
            {'multacs', 'tdma', 'ueTxPow', 10}, ...
            {'multacs', 'tdma', 'ueTxPow', 20}, ...
            {'multacs', 'fdma', 'ueTxPow', 10}, ...
            {'multacs', 'fdma', 'ueTxPow', 20} };
        legFmtStr = 'MA=%s txPow=%d dBm'; % Format string for legend
        
        
    case 'sectorization'
        % Sectorization test
        
        paramTest = { ...
            {'nantBS', 64, 'nantUE', 64, 'sectPico', 0}, ...
            {'nantBS', 64, 'nantUE', 64, 'sectPico', 1} };
        legFmtStr = 'nant=%d x %d sect=%d'; % Format string for legend
        
    case 'user'
        % # user test
        paramTest = { ...
            {'nuePerPico', 5, 'multacs', 'tdma'}, ...
            {'nuePerPico', 10, 'multacs', 'tdma'}, ...
            {'nuePerPico', 20, 'multacs', 'tdma'}, ...
            {'nuePerPico', 5, 'multacs', 'fdma'}, ...
            {'nuePerPico', 10, 'multacs', 'fdma'}, ...
            {'nuePerPico', 20, 'multacs', 'fdma'} };
        legFmtStr = 'nue/cell=%d MA=%s';
        
    case 'TxPow'
        % Spectral efficiency vs Tx pow
        
        plotAvgSE = true;
        xval = (20:5:40)';
        %xval = (10:5:30)';
        paramTest = cell(1,length(xval));
        for k = 1 : length(xval)
            paramTest{1,k} = {'picoTxPow', xval(k), 'ueTxPow', xval(k)-10};
            %paramTest{1,k} = {'picoTxPow', xval(k), 'ueTxPow', 20};
            %paramTest{1,k} = {'picoTxPow', 30, 'ueTxPow', xval(k)};
        end
        legFmtStr = 'BSTxPow = %d dBm, UETxPow = %d dBm';
        xLblStr = 'Tx Power (dBm)';
        
    case 'freq'
        % Spectral efficiency vs carrier freq using new data
        paramTest = {...
            {'freq', 28, 'nantUE', 16}, ...
            {'freq', 28, 'nantUE', 64}, ...
            {'freq', 73, 'nantUE', 16}, ...
            {'freq', 73, 'nantUE', 64} };
        legFmtStr = '%d GHz, nantUE = %d';
        
    case 'block'
        % Spectral efficiency with and without blocking
        paramTest = {...
            {'hybridDist', 'curvefit', 'blockDist', 0} ...
            {'hybridDist', 'curvefit', 'blockDist', -50}, ...
            {'hybridDist', 'curvefit', 'blockDist', -75}, ...
            {'hybridDist', 'curvefitNLOS', 'blockDist', -50}
            };
        legFmtStr = 'mod=%s, T=%d';
        
    case 'BW'
        % Spectral efficiency vs. BW
        plotAvgSE = true;
        xval = (100:300:1000)';
        paramTest = cell(1,length(xval));
        for k = 1 : length(xval)
            paramTest{1,k} = {'bwMHz', xval(k)};
        end
        legFmtStr = 'BW = %d MHz';
        xLblStr = 'Bandwidth (MHz)';
        
    case 'ISD'
        % Spectral efficiency vs. ISD
        
        plotAvgSE = true;
        defArea = 4000000; % kmsq - TODO: find a programmatic way to get these
        defUePerPico = 10;
        defPicoRad = 100;
        defNUe = defUePerPico*round(defArea/(pi*defPicoRad^2));
        
        xval = 2*(100:10:200)';     % intersite distance
        tempNbs = round(defArea./(pi*(xval/2).^2));
        nuePerPicos = defNUe./tempNbs;
        paramTest = cell(1,length(xval));
        for k = 1 : length(xval)
            paramTest{1,k} = {'picoRadius', xval(k)/2, 'nuePerPico', nuePerPicos(k)};
        end
        legFmtStr = 'picoRadius (ISD/2) = %d m, nuePerPico = %6.4f';
        xLblStr = 'Inter-site Distance (m)';
        
    case 'sect'
        paramTest = { ...
            {'picoHex', 0}, ...
            {'picoHex', 1} ...
            };
        legFmtStr = 'hex = %d';
        
    case 'LOS'
        % Power test
        paramTest = { ...
            {'plModType', 'dist'}, ...
            {'plModType', 'hybrid'} ...
            };
        legFmtStr = 'mod=%s';
        
    case 'nstream'
        % Tests effect of single stream processing at BS
        paramTest = { ...
            {'nstreamBS', 0}, ...
            {'nstreamBS', 1}, ...
            {'nstreamBS', 2}, ...
            {'nstreamBS', 4}, ...
            };
        legFmtStr = 'nstream=%d';
end

ntest = length(paramTest);



rateTestUL = cell(ntest,1);
rateTestDL = cell(ntest,1);
sinrTestUL = cell(ntest,1);
sinrTestDL = cell(ntest,1);
inrTestUL = cell(ntest,1);
inrTestDL = cell(ntest,1);
seTestUL = zeros(ntest,1);
seTestDL = zeros(ntest,1);
rateTable = zeros(ntest,4);
% Loop over test
legStr = cell(ntest,1);
for itest = 1:ntest
    
    % Set parameters
    param = paramTest{itest};
    
    % Legend string
    if (length(param) == 2)
        legStr{itest} = sprintf(legFmtStr, param{2});
    elseif (length(param) == 4)
        legStr{itest} = sprintf(legFmtStr, param{2}, param{4});
    elseif (length(param) >= 6)
        legStr{itest} = sprintf(legFmtStr, param{2}, param{4}, param{6});
    end
    
    % Run simulation
    capSim;
    
    % Save results
    if (calcUL)
        rateTestUL{itest} = rateTotUL;
        sinrTestUL{itest} = sort(sinrULs(:));
        inrTestUL{itest} = sort(inrULs(:));
        rateTable(itest,2) = avgRateUL;
        rateTable(itest,4) = edgeRateUL;
        fprintf(1,'%d %s UL:  mean=%10.4e cell-edge=%10.4e \n', ...
            itest, legStr{itest}, avgRateUL , edgeRateUL);
        seTestUL(itest) = seUL;
    end
    if (calcDL)
        rateTestDL{itest} = rateTotDL;
        sinrTestDL{itest} = sort(sinrDLs(:));
        inrTestDL{itest} = sort(inrDLs(:));
        rateTable(itest,1) = avgRateDL;
        rateTable(itest,3) = edgeRateDL;
        fprintf(1,'%d %s DL:  mean=%10.4e cell-edge=%10.4e \n', ...
            itest, legStr{itest}, avgRateDL,edgeRateDL );
        seTestDL(itest) = seDL;
    end
end
for itest = 1 : ntest
    disp([legStr{itest} '	' num2str(rateTable(itest,:))]);
end

% Plot UL rate CDF
if (ntest == 1)
    colors = 1;
else
    colors = (0:1/(ntest-1):1)';
end
if (calcUL)
    
    if (calcDL)
        subplot(1,2,1);
    else
        subplot(1,1,1);
    end
    for itest = 1 : ntest
        nue = length(rateTestUL{itest});
        p = (1:nue)/nue;
        semilogx(rateTestUL{itest}, p, '-','LineWidth',2, 'Color', [colors(itest) 1-colors(itest) 1]);
        hold on;
    end
    hold off;
    grid on;
    set(gca,'FontSize',16);
    xlabel('Rate (Mbps)');
    ylabel('Cummulative prob');
    title('Uplink')
    axis([1 1e3 0 1]);
    legend(legStr,'Location','NorthWest');
end

% Plot DL rate CDF
if (calcDL)
    
    if (calcUL)
        subplot(1,2,2);
    else
        subplot(1,1,1);
    end
    for itest = 1 : ntest
        nue = length(rateTestDL{itest});
        p = (1:nue)/nue;
        semilogx(rateTestDL{itest}, p, '-','LineWidth',2, 'Color', [colors(itest) 1-colors(itest) 1]);
        hold on;
    end
    hold off;
    grid on;
    set(gca,'FontSize',16);
    xlabel('Rate (Mbps)');
    ylabel('Cummulative prob');
    title('Downlink');
    axis([1 1e3 0 1]);
    legend(legStr);
end

if(plotAvgSE)
    figure;
    
    seLegStr ={};
    subplot(1,2,1)
    if (calcDL)
        plot(xval,seTestDL,'b-s','LineWidth',2);
        hold on;
        seLegStr = [seLegStr, 'downlink'];
    end
    if (calcUL)
        semilogy(xval,seTestUL,'g--o','LineWidth',2);
        seLegStr = [seLegStr, 'uplink'];
    end
    grid on;
    set(gca,'FontSize',16);
    xlabel(xLblStr);
    ylabel('Spectral Efficiency (bps/Hz/cell)');
    if ~isempty(seLegStr)
        legend(seLegStr,'Location','best');
    end
    subplot(1,2,2)
    if (calcDL)
        plot(xval,rateEdgeDL,'b-s','LineWidth',2);
        hold on;
    end
    if (calcUL)
        plot(xval,rateEdgeUL,'g--o','LineWidth',2);
    end
    grid on;
    set(gca,'FontSize',16);
    xlabel(xLblStr);
    ylabel('Cell-edge User Spec. Eff. (bps/Hz/user)');
    if ~isempty(seLegStr)
        legend(seLegStr,'Location','best');
    end
end

if (savedat)
    cmd = ['save data/', testName, fileSuffix, ' testName legStr paramTest rateTestDL ', ...
        'rateTestUL sinrTestDL sinrTestUL inrTestDL inrTestUL seTestUL seTestDL rateTable'];
    eval(cmd);
end

clear param param0; % not to interfere with individual run of CapSim

if (testName == tests{4})
    UL_Sinr_SD
end