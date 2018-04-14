function [ schedParam ] = schedule( nbs, Icell, sinrObj, traffic, schedOpt)
%SCHEDULE: Round robin scheduler
%   FDMA/TDM/SDMA based scheduling of DL links in a picocellular multi-BS
%   environment.
% Input:
% nbs     : number of base stations
% Icell   : UE's assicaition with BS
% se      : spectral efficienncy, should consider power splitting due to SDMA
% schedOpt: options required for scheduling
% Output:
% schedParam : Results obtained from scheduled data Tx, like rate, etc.
%              Extendable to contain different metrics.



multiacs = schedOpt.multiacs; % FDMA / SDMA, SDMA with K=1 => TDMA
K = schedOpt.numstreams;      % Number of streams at the BS, K = nant for digital BF, K = 1 for analog/tdma
numsf = schedOpt.nsf;         % total number of subframes being simulated
bw = schedOpt.bw;             % total bandwidth per BS (i.e channel BW) in MHz
eta = schedOpt.eta;           % percentage of SF not used for control/channel estimation/reverse direction 

tti = (traffic.getTTI())*1e3;        % subframe length in ms; in trafficGen 'tti' is stored in seconds

% for RR scheduler
nue = length(Icell);
bsUserIndex = ones(1,nbs);
totDoF = eta*bw*tti*1e3; % total available degrees of freedom; (MHz*ms = 1e3),per stream
% do not consider a RB based frame; have a constant Control/CSI overhead.
totDataTxue = zeros(1,nue);

% For each subframe and each base station
for sfnum = 1:numsf
    
    %usrQueue = generateTraffic(usrQueue, tti, trafficParams); % need to implement
    traffic.genTraffic(sfnum);
    
    for bs=1:nbs
        ueind = find(Icell == bs); % indices of UEs attached to bs
        %dataQueue = usrQueue(Icell == bs);
        
        if (isempty(ueind))
            continue; % no ue associated with BS
        end
        
        nuebsi = length(ueind);
        nxtIndx = bsUserIndex(bs);

        %startIndx = nxtIndx;
        
        pktSz = traffic.getQueue(ueind);
        
        if (sum(pktSz~=0) == 0)
            continue;
        end
        
        % FDMA can be done with digital BF and multiple streams
        if (strcmpi(multiacs,'fdma'))

            % Equally split the BW between all users with Data
            % Not ideal, but this is RR scheduling in OFDMA terms.
            % (Verify this)
            Dlopt.ueSched = ueind(pktSz~=0);
            Dlopt.bwSched = ones(length(Dlopt.ueSched),1)*(bw/length(Dlopt.ueSched)); % equal BW distribution
            Dlopt.powSplit = 1; % No split in total power for FDMA

            specEff = sinrObj.DlSinrCalcBSi(Dlopt);
            
            pktSz = pktSz(pktSz~=0);
            txData = min( pktSz, (eta*(specEff.*Dlopt.bwSched*1e6)*tti*1e-3)');
            
            traffic.dequeue(txData, Dlopt.ueSched);
            
            totDataTxue(Dlopt.ueSched) = totDataTxue(Dlopt.ueSched) + txData;
            
        elseif (strcmpi(multiacs,'sdma') || strcmpi(multiacs,'tdma') )
            
            Dlopt.powSplit = min(K, sum(pktSz~=0));           % Max number of beams in which the power is split

            
            
            if (Dlopt.powSplit < K)
                % there are exactly n < K users with data in the cell
                Dlopt.ueSched = ueind(pktSz~=0); % all the users are scheduled
                pktSz = pktSz(pktSz~=0);
            else
                icnt = 1;
                pkt = zeros(1,K);
                while icnt <= K
                    if (pktSz(nxtIndx) ~= 0)
                        Dlopt.ueSched(icnt) = ueind(nxtIndx);
                        pkt(icnt) = pktSz(nxtIndx);
                        icnt = icnt + 1;
                    end
                    nxtIndx = (nxtIndx == nuebsi)*1 + (nxtIndx < nuebsi)*(nxtIndx+1);

                end
                pktSz = pkt;
            end
            
            %Dlopt.ueSched = uetosched(schedIndex);
            Dlopt.bwSched = ones(length(Dlopt.ueSched),1)*bw; % assign total bandwidth

            specEff = sinrObj.DlSinrCalcBSi(Dlopt);
            
            
            
            txData = min( pktSz, (eta*(specEff.*Dlopt.bwSched*1e6)*tti*1e-3)');
            
            traffic.dequeue(txData, Dlopt.ueSched);
            
            totDataTxue(Dlopt.ueSched) = totDataTxue(Dlopt.ueSched) + txData;
            
            
%             numStreams = Dlopt.powSplit; % upper bounded by the number of users
% 
%             % schedule each stream
%             for ns=1:numStreams
%                 
%                 uei = Dlopt.ueSched(ns);
%                 
%                 pktSzi = pktSz(ns);
%                 
%                 if (pktSzi == 0)
%                     error('Programming error; debug code');
%                     % should not come here
%                     % end the loop without increasing ns
%                 end
%                 seij = specEff(ns);
%                 
%                 txData = min( eta*(seij*Dlopt.bwSched(ns)*1e6)*tti*1e-3, pktSzi);
%                 traffic.dequeue(txData, uei);
%                 totDataTxue(uei) = totDataTxue(uei) + txData;
% 
%             end
%             
%             lstServed = find(ueind == uetosched(schedIndex(end)));
%             nxtIndx = (lstServed == nuebsi)*1 + (lstServed < nuebsi)*(lstServed+1);

        else
            error('Unknown/unsupported multiple access scheme');
        end
        
        bsUserIndex(bs) = nxtIndx;

   end
end

totDataQueue = traffic.getTotDataQueues();

[Eq, waitTime] = traffic.calcAvgQueueLength(numsf);

totDataTxue(totDataQueue == 0) = []; % do not consider these; no data were Txed for these UEs
Eq(totDataQueue == 0) = [];
waitTime(totDataQueue == 0) = [];
totDataQueue(totDataQueue == 0) = [];


schedParam.rate = totDataTxue/(numsf*tti*1e3);  % in Mbps
schedParam.service = totDataTxue./totDataQueue; % [0,1] --> how much was served

schedParam.avgQueueLen = Eq;
schedParam.avgWait = waitTime;



% we should do something with queue length etc.

end
    



