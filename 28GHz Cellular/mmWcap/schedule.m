function [ schedParam ] = schedule( nbs, Icell, se, traffic, schedOpt)
%SCHEDULE Summary of this function goes here
%   Detailed explanation goes here
% Input:
% nbs     : number of basestations
% Icell   : UE's assicaition with BS
% se      : spectral efficienncy, should consider power splitting due to SDMA
% schedOpt: options required for scheduling



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
        ueind = find(Icell == bs);
        %dataQueue = usrQueue(Icell == bs);
        
        if (isempty(ueind))
            continue; % no ue associated with BS
        end
        
        nuebsi = length(ueind);
        
        nxtIndx = bsUserIndex(bs);
        %availDof = totDoF;
        startIndx = nxtIndx;
        % FDMA can be done with digital BF and multiple streams
        if (strcmpi(multiacs,'fdma'))
            availDof = totDoF; % check this. Should this be DoF*K? Discuss.
            
            while (availDof > 0)
                uei = ueind(nxtIndx);
                seij = se(uei);
                pktSz = traffic.getQueue(uei);
                
                dofReq = pktSz/seij;
                
                if (dofReq > availDof)
                    txData = availDof*seij;
                    availDof = 0;
                else
                    txData = pktSz;
                    availDof = availDof - dofReq;
                end
                traffic.dequeue(txData, uei);
                
                totDataTxue(uei) = totDataTxue(uei) + txData;
                %Update the next UE index to be scheduled
                nxtIndx = (nxtIndx == nuebsi)*1 + (nxtIndx < nuebsi)*(nxtIndx+1);
                if (startIndx == nxtIndx) 
                    break; % exhausted all set of users.
                end
                
            end
            
        elseif (strcmpi(multiacs,'sdma') || strcmpi(multiacs,'tdma') )
            availDof = totDoF;
            
            % schedule each stream
            for schedInd = 1:K
                uei = ueind(nxtIndx);
                seij = se(uei);
                pktSz = traffic.getQueue(uei);
                
                dofReq = pktSz/seij;
                if (dofReq > availDof)
                    txData = availDof*seij;
                else
                    txData = pktSz;
                end
                
                traffic.dequeue(txData, uei);
                totDataTxue(uei) = totDataTxue(uei) + txData;
                %Update the next UE index to be scheduled
                nxtIndx = (nxtIndx == nuebsi)*1 + (nxtIndx < nuebsi)*(nxtIndx+1);
                if (startIndx == nxtIndx) 
                    break; % exhausted all set of users. When K > Nue
                end
            end
        else
            error('Unknown/unsupported multiple access scheme');
        end
        
        bsUserIndex(bs) = nxtIndx;

   end
end

totDataQueue = traffic.getTotDataQueues();
totDataTxue(totDataQueue == 0) = []; % do not consider these; no data were Txed for these UEs
totDataQueue(totDataQueue == 0) = [];
schedParam.rate = totDataTxue/(numsf*tti*1e3);  % in Mbps
schedParam.service = totDataTxue./totDataQueue; % [0,1] --> how much was served

% we should do something with queue length etc.

end
    



