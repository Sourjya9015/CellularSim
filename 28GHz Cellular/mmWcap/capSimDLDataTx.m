% Uplink and Downlink capacity simulation

% Set path to contain all the subdirectories of NetSim
addpath(genpath('..'));

% Parameters
picoRadius = 100;  % cell radius for picocells
nuePerPico = 20;   % num ues per macro sector
posLim0 = [1000 1000]';     % area for region in meters
chanIfName = 'lte-chan-if'; % name for all the channel interfaces
devName = 'lte-dev';        % name for all the LTE network devices
picoTxPow = 35;     % Pico TX power in dBm
ueTxPow = 20;       % UE TX power in dBm
picoNoiseFig = 5;   % Pico noise figure
ueNoiseFig = 10;     % UE Total (accross all Nrx antennas)noise figure
nantBS = 64;        % num antennas at the BS
nantUE = 16;        % num antennas at the UE
fcMHz = 28e3;       % carrier freq in MHz
seMax = 4.8;        % max spec. efficiency
sinrLossdB = 3;     % loss relative to Shannon capacity
bwMHz = 1000;       % Total bandwidth in MHz
picoHex = true;     % Picocells sectorized with hex layout
nsectPico = 3;      % num pico sectors in hex layout
multacs = 'tdma';   % FDMA / TDMA flag
multacsDL = 'sdma'; % TDMA/FDMA/SDMA
bsNumStreams = 4;   % num of streams; equal to the num of RF streams
ndrop = 3;          % number of drops
calcUL = false;     % calculate UL capacity
calcDL = true;      % calculate DL capacity
plModType = 'hybrid';    % select either 'hybrid' or 'dist'
% if 'hybrid', then set parameters below
hybridDist = 'curvefitNLOS'; % hybrid distribution type (see below)
blockDist = 0;      % blocking distance
nstreamBS = 0;      % max num streams at BS (0 = infinite)
freq = 28;          %GHz %%%% 0 = use old default PL model
speedoflight = 0.3; % Gm/s
%pl0LOS = 81.4;      % path loss at d=1m for LOS model
% leave empty to calculate wrto Friis' Law
intNull = false;    % enable interference nulling

pedge = 0.05;       % cell-edge percentage = 5
dutyCycle = 0.5;    % 50 percent is for DL
overhead = 0.2;     % 20 percent overhead

nbits = 4;          % quantization bits
[~,alpha,~] = unifQuant(nbits);

% Process parameters for batch if a parameter string was supplied
if exist('param0', 'var')
    nparam = length(param0)/2;
    for iparam = 1:nparam
        cmd = sprintf('%s = param0{2*iparam};', param0{2*iparam-1});
        eval(cmd);
    end
end
if exist('param', 'var')
    nparam = length(param)/2;
    for iparam = 1:nparam
        cmd = sprintf('%s = param{2*iparam};', param{2*iparam-1});
        eval(cmd);
    end
end

% Get NLOS pathloss model based on freq
if (freq == 28)
    %load ../chanMod/data/chanParam28.mat
    a1 = 0.0334;
    b1 = 5.2;
    a2 = 0.0149;
    plexpNLOS = 2.92;
    pl0NLOS = 72;
    shadStdNLOS = 8.7;
    
    plexpLOS = 2;       % path loss explonent for LOS %TODO: update for 28 
    shadStdLOS = 6.0;     
elseif (freq == 73)
    %load ../chanMod/data/chanParam73.mat
    a1 = 0.0334;%0.0305;
    b1 = 5.2;%4.55;
    a2 = 0.0149;%0.028;
    plexpNLOS = 2.69;
    pl0NLOS = 82.7;
    shadStdNLOS = 7.7;
    
    plexpLOS = 2;       % path loss explonent for LOS
    shadStdLOS = 5.8;
else
    pl0NLOS = 75.85 - 8.56; % path loss at d=1m for total path
    plexpNLOS = 3.73;       % path loss exponent for NLOS
    shadStdNLOS = 5.56;       % per path shadowing
    a1 = 0.0359;
    b1 = 5.57;
    a2 = 0.0210;
end
% if (freq ~= 0)
%     pl0NLOS = plcoeff(1);
%     plexpNLOS = plcoeff(2);
%     shadStdNLOS = 5.56;
% end
pl0LOS = 20*log10(4*pi*freq/speedoflight);

% Set hybrid probability model
plModOpt.plModType = plModType;
switch hybridDist
    case 'NLOS'
        probR = @(R)zeros(size(R));   % only NLOS
        poutR = @(R)zeros(size(R));
    case 'LOS'
        probR = @(R)ones(size(R));    % only LOS
        poutR = @(R)zeros(size(R));
    case 'relay-ue'                   % 3GPP case1 relay-UE
        probR =@(R)0.5-min(0.5,5*exp(-156./R))+min(0.5,5*exp(-R/30));
        poutR = @(R)zeros(size(R));
    case 'bs-ue'                      % 3GPP case1 BS-UE
        probR =@(R)min(18./R,1).*(1-exp(-R/63))+exp(-R/63);
        poutR = @(R)zeros(size(R));
    case 'block'
        thetaBlk = 2*0.256;
        poutR = @(R)1-(min(18./(thetaBlk*R),1).*(1-exp(-(thetaBlk*R)/63))+exp(-(thetaBlk*R)/63));
        %pl0LOS = 200;  not needed anymore % Very high number to block the transmissions
        probR = @(R)zeros(size(R)); % all NLOS
    case 'blockT'
        poutR = @(R)(R>blockDist);
        % pl0LOS = 200; not needed anymore % Very high number to block the transmissions
        probR = @(R)zeros(size(R)); % all NLOS
    case 'curvefit'
        probR = @(R)exp(-a2*R);
        poutR = @(R)(1-min(1,exp(-a1*(R-blockDist)+b1)));
    case 'curvefitNLOS'
        probR = @(R)zeros(size(R));
        poutR = @(R)(1-min(1,exp(-a1*(R-blockDist)+b1)));
    case 'curvefitLOS'
        probR = @(R)ones(size(R));    % only LOS
        poutR = @(R)(1-min(1,exp(-a1*(R-blockDist)+b1)));
end

% Determine layout of pcios.  Also recompute the exact position limit to
% insure that they agree with the position limits
if (picoHex)
    picoLM = HexLayout(2*picoRadius, posLim0);
    [nx,ny] = picoLM.getDim();
    npicoNode = nx*ny;
    npico  = npicoNode*nsectPico;
    posLim = picoLM.getPosLim();
else
    npicoNode = round( posLim0(1)*posLim0(2)/(pi*picoRadius^2) );
    picoLM = UnifLayout(posLim0);
    posLim = posLim0;
    npico = npicoNode;
end

% Compute number of other elements
nue = round( npico*nuePerPico );

% BF gain computation options
nantBS1 = round(sqrt(nantBS));
nantUE1 = round(sqrt(nantUE));
bfOpt.intNull = intNull;
bfOpt.covFn = sprintf('../chanMod/data/covData%d_%dx%dx%dx%d', ...
    freq, nantBS1, nantBS1, nantUE1, nantUE1);
%nb
bfOpt.nantBS = nantBS;

% Create node data
picoDat.num = npicoNode;
ueDat.num   = nue;
picoDat.layoutMod   = picoLM;
ueDat.layoutMod     = UnifLayout(posLim);
nodeDatMap = containers.Map();
nodeDatMap('pico') = picoDat;
nodeDatMap('ue') = ueDat;

% Create channel data
chanDat = [];
chanName = 'mmW';
chanDatMap = containers.Map();
chanDatMap(chanName) = chanDat;

% Create channel interface data
% In this case, we have one interface for each node type
nodeNames = nodeDatMap.keys();
chanIfDatMap = containers.Map();
nnode = length(nodeNames);
for it = 1:nnode
    nodeName = nodeNames{it};
    chanIfDat.chanName = chanName;
    chanIfDat.nodeName = nodeName;
    if (strcmp('pico',nodeName) && picoHex)
        chanIfDat.nsect = nsectPico;
    else
        chanIfDat.nsect = 1;
    end
    ifName = strcat(nodeName,':',chanName);
    chanIfDatMap(ifName) = chanIfDat;
end

% This clear is needed due to a MATLAB bug that calls a delete on chan
% when there is still an outstanding reference in a map object.
clear net chan;

% Create the het net
net = HetNet(nodeDatMap, chanDatMap, chanIfDatMap);

% Set global network properties
net.set('posLim', posLim);

% Set channel parameters
chan = net.getChan(chanName);
chan.set('bwMHz', bwMHz, 'fcMHz', fcMHz);

% Set the channel interface parameters
picoChanIf  = net.getChanIfList('pico:mmW');
ueChanIf    = net.getChanIfList('ue:mmW');

% Set device specific interface parameters
for ipico = 1:npico
    chanIf = picoChanIf{ipico};
    chanIf.set('txPowdBm', picoTxPow, 'noiseFig', picoNoiseFig);
    if (picoHex)
        ap = SectAntPattern();
        ap.set('patType', 'flat', 'sectInd', chanIf.sectInd,'gainMax',0,...
            'fbGain',80,'secWid',2*pi/nsectPico);
        chanIf.addAntPattern(ap);
    end
end
for iue = 1:nue
    ueChanIf{iue}.set('txPowdBm', ueTxPow, 'noiseFig', ueNoiseFig);
end

% Add path loss models.
net.addPlMod('picoPL',  'pico:mmW',  'ue:mmW', plModOpt);

% Set the path loss params
picoPL = net.getPlMod('picoPL');
switch plModOpt.plModType
    case 'hybrid'
        picoPL.set('plexpLOS', plexpLOS, 'plexpNLOS', plexpNLOS, ...
            'shadStdLOS', shadStdLOS, 'shadStdNLOS', shadStdNLOS,...
            'pl0LOS', pl0LOS, 'pl0NLOS', pl0NLOS,...
            'probR', probR, 'poutR', poutR);
        
    case 'dist'
        picoPL.set('plexp', plexpNLOS,'shadStd', shadStdNLOS, 'pl0', pl0NLOS);
end

% Loop over drops
if (calcUL)
    rateULs = zeros(nue,ndrop);
    sinrULs = zeros(nue,ndrop);
    inrULs = zeros(nue,ndrop);
    seUL = 0;
end
if (calcDL)
    rateDLs = zeros(nue,ndrop);
    sinrDLs = zeros(nue,ndrop);
    inrDLs = zeros(nue,ndrop);
    seDL = 0;
    %nb 
    sinrDLeffs = zeros(nue,ndrop);
    Ps = zeros(nue,ndrop);
    N0 = zeros(nue,ndrop);
    bfGainRxDLs = zeros(nue,ndrop);
    Pintfr = zeros(nue,ndrop);
    bwmq = zeros(nue,ndrop);
    
    schedRateDL = [];
    serviceDL = [];
end

tti = 1e-3; % in seconds
nsf = 1000;
trafficTime = tti*nsf;
% Set traffic parameters
trafficOpt.lambda = [100  10  1] ; % packets/sec
trafficOpt.size = [1e3 2e5 1e9];  % mean size in Bytes
trafficOpt.frac = [0.3 0.3 0.4];  % 80% users in small packets, 20% large packets
trafficOpt.nqueue = nue;
trafficOpt.tti = tti; % 1 ms, check 3GPP spec for this
trafficOpt.totTime = trafficTime; % in second(s)

trafficHandle = trafficGen(trafficOpt);

for idrop = 1 : ndrop
    fprintf(1,'Drop %d of %d\n', idrop, ndrop);
    % Random drop
    net.drop();
    %for i=1:ntti
    % Get path loss and DL SINR
    [sinrDL, Icell, pathLoss, txpowDL] = chan.computeSinr('ue:mmW', {'pico:mmW'}); % args: rx,tx
    % Notes SD: 
    % pathLoss -> 390x3900 matrix tabling the pathloss between each UE-BS pair. 
    % Icell    -> 1x3900   vector recording BS each UE is associated to. 
    % txpowDL  -> 1x390    DL Tx power of BSs.
    
    % Apply the BF gain -- assuming all links have a spatial structure as
    % the NLOS channel
    bfOpt.txpow = repmat( picoTxPow, npico, 1);    % max TX pow in dBm
    bfOpt.noisepow = chan.kT + 60 + 10*log10(bwMHz) + ueNoiseFig;
    [pathLossDL,pathLossUL,bfGainRxDLdes] = applyBfGain(pathLoss', Icell, bfOpt);
    display('Checking size: bfGainRxDLdes:');
    size(bfGainRxDLdes)
    % Reduce BF gain if BS is reduced to single stream processing
    if (nstreamBS > 0)
        bfOpt.nstream = nstreamBS;
        [pathLossDL,pathLossUL] = bfGainRedStream( pathLossDL, pathLossUL, Icell, bfOpt);
    end
    
    % Uplink capacity estimation
    % ----------------------------
    if (calcUL)
        opt.pathLoss = pathLossUL;              % path loss in dB
        opt.txpow = repmat( ueTxPow, nue, 1);   % max TX pow in dBm
        opt.Icell = Icell;                      % Icell(j) = RX index for TX i
        if strcmp(multacs,'fdma')
            opt.fdma = true;                    % Is FDMA is enabled?  Otherwise, assume TDMA
        else
            opt.fdma = false;
        end
        opt.bwMHzTot = bwMHz;   % Total bandwidth in MHz
        opt.noisepow = chan.kT + 60 + 10*log10(bwMHz) + picoNoiseFig;
        % noise power in dBm
        
        % Find optimal TX power
        ulOptFn = ULPowConOptFn(opt);
        p = ulOptFn.optPower();
        
        % Get rate and SINR distributions
        [rateUL,sinrUL,inrUL] = ulOptFn.computeRate(p);
        mySNIRUL = sinrUL;
        rateULs(:,idrop) = rateUL;
        sinrULs(:,idrop) = 10*log10(sinrUL);
        inrULs(:,idrop) = 10*log10(inrUL);
        seUL = seUL+ sum(rateUL)/npico/bwMHz;
        
    end
    
    % Downlink capacity estimation
    % ----------------------------
    if (calcDL)
        opt.pathLoss = pathLossDL;                                          % path loss in dB
        opt.txpow = repmat( picoTxPow, npico, 1);                           % max TX pow in dBm
        opt.Icell = Icell;                                                  % Icell(j) = RX index for TX i
        opt.bwMHzTot = bwMHz;                                               % Total bandwidth in MHz
        opt.noisepow = chan.kT + 60 + 10*log10(bwMHz) + ueNoiseFig;         % for M antennas
        opt.rxBFgains = bfGainRxDLdes;                                      % Rx BF gains to the home BS
        opt.alpha = alpha;              

        
        % nultiple access based parameters
        multiacsOpt.multiacs = multacsDL;
        multiacsOpt.K = bsNumStreams;
        
        % Calculate rates
        dlSinrCalc = DLSinrCalc(opt, multiacsOpt);
        rateDL = dlSinrCalc.rate;     % all these are per UE (3900x1)
        sinrDL = dlSinrCalc.sinr;   %nb can change here or inside!
        % nb: sinrDL = 
        
        inrDL = dlSinrCalc.inr;
        seDL = dlSinrCalc.specEff;
        
        % Do scheduler/traffic generator
        % Add quantized SNR
        schedOpt.multiacs = multacsDL;
        schedOpt.numstreams = bsNumStreams;
        schedOpt.nsf = nsf;   % number fo SFs
        schedOpt.bw  = bwMHz; % in MHz
        schedOpt.eta = dutyCycle*(1-overhead);
        
        schedOut = schedule( npico, Icell, dlSinrCalc, trafficHandle, schedOpt);
        
        schedRateDL = [schedRateDL, schedOut.rate];
        serviceDL = [serviceDL, schedOut.service];
        
        %nb
        Ps(:,idrop) = dlSinrCalc.psign;
        N0(:,idrop) = dlSinrCalc.pnoise;
        Pintfr(:,idrop) = dlSinrCalc.pintrfr;
        
        rateDLs(:,idrop) = rateDL;
        sinrDLs(:,idrop) = sinrDL;
        inrDLs(:,idrop) = inrDL;
        seDL = seDL + sum(rateDL)/npico/bwMHz;
        bfGainRxDLs(:,idrop) = bfGainRxDLdes;
        
        %nb        
        bwmq(:,idrop) = dlSinrCalc.bwMHz;        
 
        
        
    end
end
if (calcUL)
    seUL = seUL/ndrop;
    rateTotUL = sort(rateULs(:));
    avgRateUL = sum(rateTotUL)/(ndrop*npico)*(1-dutyCycle)*(1-overhead);
    edgeRateUL = rateTotUL(round(pedge*length(rateTotUL)))*(1-dutyCycle)*(1-overhead);
end
if (calcDL)
    seDL = seDL/ndrop;
    rateTotDL = sort(rateDLs(:));
    avgRateDL = sum(rateTotDL)/(ndrop*npico)*dutyCycle*(1-overhead);
    edgeRateDL = rateTotDL(round(pedge*length(rateTotDL)))*dutyCycle*(1-overhead);
    sinrDL = 10.^(0.1*sinrDLs(:));
    
    schedRateDL = sort(schedRateDL(:)); % per user rate
    serviceDL = sort(serviceDL(:));     % service provided per UE
    %nb
    bwmq = bwmq(:);
    snrDL = Ps./N0;
    snrDL = snrDL(:);
    
    switch hybridDist
        case 'curvefit'
            save('data/sinrDL.mat','sinrDL');
        case 'curvefitLOS'
            save('data/sinrDLlos.mat','sinrDL');
        case 'curvefitNLOS'
            save('data/sinrDLnlos.mat','sinrDL');
    end
end

% Plot results if the param cell array was not specified
% If the param cell array was specified the results will be plotted after
% running the batch test capTest
if ~exist('param','var')
    
    
    if (calcDL)
        if (calcDL && calcUL)
            subplot(1,2,1);
        else
            subplot(1,1,1);
        end
        %rateTotDL = sort(rateDLs(:));
        figure (1);
        n = length(rateTotDL);
        p = (1:n)/n;
        h = semilogx(sort(rateTotDL),p,'-');
        grid on;
        set(h,'LineWidth',2);
        set(gca,'FontSize',16);
        xlabel('Rate (Mbps)');
        ylabel('Cummulative prob');
        grid on; hold on;
        axis([1 1e4 0 1]);
        fprintf(1,'DL:  mean=%10.4e cell-edge=%10.4e \n', avgRateDL, edgeRateDL );
        
        
        % for scheduled rate
        nsched = length(serviceDL);
        prob = (1:nsched)/nsched;
        figure(2); 
        hold on;
        
        cdf = serviceDL;
        
        plot(cdf,prob,'-','Linewidth',2);
        grid on;
        axis([0 1 0 1]);
        set(gca,'FontSize',16);
        xlabel(' fraction of data delevered');
        ylabel('CDF');
        title('Service with RR Scheduler');
        %legend('Theoretical','Scheduled');
        
        
    end
    
    if (calcUL)
        if (calcDL && calcUL)
            subplot(1,2,2);
        else
            subplot(1,1,1);
        end
        %rateTotUL = sort(rateULs(:));
        n = length(rateTotUL);
        p = (1:n)/n;
        figure (1);
        h = semilogx(sort(rateTotUL),p,'-');
        grid on;
        set(h,'LineWidth',2);
        set(gca,'FontSize',16);
        xlabel('Rate (Mbps)');
        ylabel('Cummulative prob');
        grid on;
        %axis([0.1 1e3 0 1]);
        fprintf(1,'UL:  mean=%10.4e cell-edge=%10.4e \n', avgRateUL, edgeRateUL );
    end
    
end




