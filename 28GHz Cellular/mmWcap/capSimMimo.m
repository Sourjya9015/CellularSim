% Uplink and Downlink capacity simulation

% Set path to contain all the subdirectories of NetSim
addpath(genpath('..'));

% Parameters
picoRadius = 100;   % cell radius for picocells
nuePerPico = 10;   % num ues per macro sector
posLim0 = [2000 2000]';     % area for region in meters
chanIfName = 'lte-chan-if'; % name for all the channel interfaces
devName = 'lte-dev';        % name for all the LTE network devices
picoTxPow = 30;     % Pico TX power in dBm
ueTxPow = 20;       % UE TX power in dBm
picoNoiseFig = 5;   % Pico noise figure
ueNoiseFig = 7;     % UE noise figure
nantBS = 64;        % num antennas at the BS
nantUE = 64;        % num antennas at the UE
fcMHz = 28e3;       % carrier freq in MHz
seMax = 4.8;        % max spec. efficiency
sinrLossdB = 3;     % loss relative to Shannon capacity
bwMHz = 1000;       % Total bandwidth in MHz
picoHex = true;     % Picocells sectorized with hex layout
nsectPico = 3;      % num pico sectors in hex layout
multacs = 'fdma';   % FDMA / TDMA flag
ndrop = 3;          % number of drops
calcUL = true;     % calculate UL capacity
calcDL = true;      % calculate DL capacity
pl0NLOS = 75.85 - 8.56; % path loss at d=dref for total path
plexpNLOS = 3.73;       % path loss exponent for NLOS
pl0LOS = 81.4;        % path loss at d=dref
                    % leave empty to calculate wrto Friis' Law
plexpLOS = 2;       % path loss explonent for LOS
plstd = 5.56;       % per path shadowing
dref = 1;           % reference distance (in m) to calculate pl0s
plModType = 'mimo';      % select either 'hybrid' or 'dist'
                         % if 'hybrid', then set parameters below
hybridDist = 'relay-ue'; % hybrid distribution type (see below)

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

% Set hybrid probability model
plModOpt.plModType = plModType;
switch hybridDist
    case 'NLOS'
        probR = @(R)zeros(size(R));   % only NLOS
    case 'LOS' 
        probR = @(R)ones(size(R));    % only LOS       
    case 'relay-ue'                   % 3GPP case1 relay-UE   
        probR =@(R)0.5-min(0.5,5*exp(-156./R))+min(0.5,5*exp(-R/30)); 
    case 'bs-ue'                      % 3GPP case1 BS-UE  
        probR =@(R)min(18./R,1).*(1-exp(-R/63))+exp(-R/63);  
end

% Determine layout of picos.  Also recompute the exact position limit to
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
nantBS1 = round(sqrt(nantBS));      % num antennas per dimension
nantUE1 = round(sqrt(nantBS));      % num antennas per dimension
for ipico = 1:npico
    chanIf = picoChanIf{ipico};
    chanIf.set('txPowdBm', picoTxPow, 'noiseFig', picoNoiseFig, ...
        'nantH', nantBS1, 'nantV', nantBS1);
    if (picoHex)        
        ap = SectAntPattern();
        ap.set('patType', 'flat', 'sectInd', chanIf.sectInd,'gainMax',0,...
            'fbGain',80,'secWid',2*pi/nsectPico);
        chanIf.addAntPattern(ap);
    end   
end
for iue = 1:nue
    ueChanIf{iue}.set('txPowdBm', ueTxPow, 'noiseFig', ueNoiseFig, ...
       'nantH', nantUE1, 'nantV', nantUE1);
end

% Add path loss models.
net.addPlMod('picoPL',  'pico:mmW',  'ue:mmW', plModOpt);

% Set the path loss params
picoPL = net.getPlMod('picoPL');
switch plModOpt.plModType
    case 'hybrid'
        picoPL.set('plexpLOS', plexpLOS, 'plexpNLOS', plexpNLOS, ...
            'shadStd', plstd, 'pl0LOS', pl0LOS, 'pl0NLOS', pl0NLOS,...
            'dref',dref, 'fcMHz', fcMHz, 'probR', probR);
        
    case 'dist'
        picoPL.set('plexp', plexpNLOS, ...
            'shadStd', plstd, 'pl0', pl0NLOS,'dref',dref, 'fcMHz', fcMHz);
    case 'mimo'
        picoPL.set('plexp', plexpNLOS, 'shadStd', plstd, 'pl0', pl0NLOS);
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
end
for idrop = 1 : ndrop
    fprintf(1,'Drop %d of %d\n', idrop, ndrop);
    
    % Random drop
    net.drop();
    
    % Get path loss and DL SINR
    [sinrDL, Isel, pathLoss, txpowDL] = chan.computeSinr('ue:mmW', {'pico:mmW'});
    
    if 0    
    switch plModOpt.plModType
        case 'dist'
            [pathLossDL,pathLossUL] = applyBfGain(pathLoss', Isel, bfOpt);
        case 'hybrid'
            % Apply path loss with BF cosidering LOS or NLOS regions
            reg = picoPL.get('region');
            [pathLossDL,pathLossUL] = applyBfGain(pathLoss', Isel, bfOpt, reg');
    end
    end
    
    % Uplink capacity estimation
    % ----------------------------
    if (calcUL)
        opt.pathLoss = pathLossUL;                 % path loss in dB
        opt.txpow = repmat( ueTxPow, nue, 1);    % max TX pow in dBm
        opt.Isel = Isel;        % Isel(j) = RX index for TX i
        if strcmp(multacs,'fdma')
            opt.fdma = true;        % Whether FDMA is enabled.  Otherwise, assume TDMA
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
        rateULs(:,idrop) = rateUL;
        sinrULs(:,idrop) = 10*log10(sinrUL);
        inrULs(:,idrop) = 10*log10(inrUL);
        seUL = seUL+ sum(rateUL)/npico/bwMHz;
    end
    
    % Downlink capacity estimation
    % ----------------------------
    if (calcDL)
        opt.pathLoss = pathLoss;                    % path loss in dB
        opt.txpow = repmat( picoTxPow, npico, 1);    % max TX pow in dBm
        opt.Isel = Isel;        % Isel(j) = RX index for TX i
        opt.bwMHzTot = bwMHz;   % Total bandwidth in MHz
        opt.noisepow = chan.kT + 60 + 10*log10(bwMHz) + ueNoiseFig;
        opt.chan = chan;
        opt.txName = 'pico:mmW';
        opt.rxName = 'ue:mmW';
        return
        % noise power in dBm
        
        % Calculate rates
        dlSinrCalc = DLSinrCalc(opt);
        rateDL = dlSinrCalc.rate;
        sinrDL = dlSinrCalc.sinr;
        inrDL = dlSinrCalc.inr;
        rateDLs(:,idrop) = rateDL;
        sinrDLs(:,idrop) = sinrDL;
        inrDLs(:,idrop) = inrDL;
        seDL = seDL + sum(rateDL)/npico/bwMHz;
    end
end
if (calcUL)
    seUL = seUL/ndrop;
end
if (calcDL)
    seDL = seDL/ndrop;
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
        rateTotDL = sort(rateDLs(:));
        n = length(rateTotDL);
        p = (1:n)/n;
        h = semilogx(sort(rateTotDL),p,'-');
        grid on;
        set(h,'LineWidth',2);
        set(gca,'FontSize',16);
        xlabel('Rate (Mbps)');
        ylabel('Cummulative prob');
        grid on;
        axis([0.1 1e3 0 1]);
    end
    
    if (calcUL)
        if (calcDL && calcUL)
            subplot(1,2,2);
        else
            subplot(1,1,1);
        end
        rateTotUL = sort(rateULs(:));
        n = length(rateTotUL);
        p = (1:n)/n;
        h = semilogx(sort(rateTotUL),p,'-');
        grid on;
        set(h,'LineWidth',2);
        set(gca,'FontSize',16);
        xlabel('Rate (Mbps)');
        ylabel('Cummulative prob');
        grid on;
        axis([0.1 1e3 0 1]);
    end
    
end

