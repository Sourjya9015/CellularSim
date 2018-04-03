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
        
        if (isempty(pktSz~=0))
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
            
            Dlopt.ueSched = ueind;
            Dlopt.bwSched = ones(length(Dlopt.ueSched),1)*bw; % assign total bandwidth
            Dlopt.powSplit = min(K, sum(pktSz~=0));           % Max number of beams in which the power is split
            
            specEff = sinrObj.DlSinrCalcBSi(Dlopt);
            
            
            numStreams = Dlopt.powSplit; % upper bounded by the number of users
            
            ns = 1;
            startIndex = nxtIndx;
            flag = 0;

            % schedule each stream
            
            while (ns <= numStreams)
                
                % Prevents infinite loops when no data is available or
                % number of streams available is very large
                if (nxtIndx == startIndex)
                    flag = flag + 1;
                end
                
                if (flag > 1) 
                    break; 
                end
                %------------------------------------------------------
                
                uei = Dlopt.ueSched(nxtIndx);
                
                pktSzi = pktSz(nxtIndx);
                
                if (pktSzi == 0)
                    nxtIndx = (nxtIndx == nuebsi)*1 + (nxtIndx < nuebsi)*(nxtIndx+1);
                    continue; % end the loop without increasing ns
                end
                seij = specEff(nxtIndx);
                
                txData = min( eta*(seij*Dlopt.bwSched(nxtIndx)*1e6)*tti*1e-3, pktSzi);
                traffic.dequeue(txData, uei);
                totDataTxue(uei) = totDataTxue(uei) + txData;
                
                nxtIndx = (nxtIndx == nuebsi)*1 + (nxtIndx < nuebsi)*(nxtIndx+1);
                ns = ns+1;
            end
%             for schedInd = 1:numStreams
%                 uei = ueind(nxtIndx);
%                 seij = specEff(nxtIndx);
%                 pkt = pktSz(nxtIndx);
%                 
%                 dofReq = pktSz/seij;
%                 if (dofReq > availDof)
%                     txData = availDof*seij;
%                 else
%                     txData = pktSz;
%                 end
%                 
%                 traffic.dequeue(txData, uei);
%                 totDataTxue(uei) = totDataTxue(uei) + txData;
%                 %Update the next UE index to be scheduled
%                 nxtIndx = (nxtIndx == nuebsi)*1 + (nxtIndx < nuebsi)*(nxtIndx+1);
%                 if (startIndx == nxtIndx) 
%                     break; % exhausted all set of users. When K > Nue
%                 end
%             end
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
    



