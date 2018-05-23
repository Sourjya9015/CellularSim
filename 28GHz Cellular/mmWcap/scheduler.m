classdef scheduler < hgsetget
    %SCHEDULER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        multiAccess = 'tdma';
        nStreams = 1;
        numSF;  % num subframes being simulated
        bwMHz;  % bandwidth
        eta;    % overhead; due to control etc.
        ttims;  % TTI in milliseconds
        
        
        trafficObj;     % An object of the traffic class
        nBS;            % Number of BS in the simulation
        Icell;          % UE associations with each BS
        DLSinrObj;      % The DL SINR calculation object
        
        thrHistory;
        
        sdmaICIReject = 1;
        
    
    end
    
    properties (Access='private')
        totDataTx;
        numUsrSched;
        fullBufferRate;
        
        sinrRecord;
        schedCount;
    end

    methods (Access='public')
        
        function obj = scheduler(opt)
            % this can be also done using set functions
            if (nargin ~=0)
                obj.nStreams = opt.nStream;
                obj.multiAccess = opt.multiAccess;
                obj.numSF = opt.numSF;
                obj.bwMHz = opt.bwMHz;
                obj.eta = opt.eta;
                obj.nBS = opt.nBS;
                %obj.Icell = opt.Icell;
                obj.ttims = opt.ttims;
            end
        end
        
        function setTrafficModel(obj,traffic)
            obj.trafficObj = traffic;
        end
        
        function setSINRModel(obj,sinrObj)
            obj.DLSinrObj = sinrObj;
        end
        
        function setCellAllocation (obj, CellSel)
            obj.Icell = CellSel;
        end
        
        function schedOpt = schedule(obj)
            
            obj.thrHistory = ones(1, length(obj.Icell))/length(obj.Icell); % Initialized to 1/Nue; 
            % Do not initialize to 0 as that will lead to possible division
            % errors!
            
            obj.totDataTx = zeros(1, length(obj.Icell));
            obj.numUsrSched = zeros(obj.nBS, obj.numSF);
            
            obj.fullBufferRate = zeros(1, length(obj.Icell));
            obj.sinrRecord = zeros(1, length(obj.Icell));
            obj.schedCount = zeros(1, length(obj.Icell));
            
            for sfnum = 1:obj.numSF
                
                obj.trafficObj.genTraffic(sfnum);
                
                for bsi = 1:obj.nBS
                    
                    ueind = find(obj.Icell == bsi); % indices of UEs attached to bs

                    if (isempty(ueind))
                        continue; % no ue associated with BS
                    end

                    %nuebsi = length(ueind);
                    %nxtIndx = bsUserIndex(bs);

                    pktSz = obj.trafficObj.getQueue(ueind);

                    if (sum(pktSz) == 0)
                        continue;
                    end
                    
                    if (strcmpi(obj.multiAccess,'tdma'))
                        obj.schedTDMA(ueind, sfnum, bsi);

                    elseif (strcmpi(obj.multiAccess,'sdma'))
                        obj.schedSDMA(ueind, sfnum, bsi);

                    elseif (strcmpi(obj.multiAccess,'fdma'))
                        obj.schedFDMA(ueind, sfnum, bsi);

                    end
            
                end
            end
            
            totDataQueue = obj.trafficObj.getTotDataQueues();

            [Eq, waitTime] = obj.trafficObj.calcAvgQueueLength(obj.numSF);

            obj.totDataTx(totDataQueue == 0) = []; % do not consider these; no data were Txed for these UEs
            Eq(totDataQueue == 0) = [];
            waitTime(totDataQueue == 0) = [];
            obj.fullBufferRate(totDataQueue == 0) = [];
            totDataQueue(totDataQueue == 0) = [];
            
            obj.sinrRecord(obj.schedCount == 0) = [];
            obj.sinrRecord = obj.sinrRecord./obj.schedCount;


            schedOpt.rate = obj.totDataTx/(obj.numSF*obj.ttims*1e3);  % in Mbps (1e-3*1e6 in the denom)
            schedOpt.service = obj.totDataTx./totDataQueue; % [0,1] --> how much was served

            schedOpt.avgQueueLen = Eq;
            schedOpt.avgWait = waitTime;
            schedOpt.fullBuffRate = obj.fullBufferRate/(obj.numSF*obj.ttims*1e3);
            schedOpt.sinr = obj.sinrRecord; % in linear scale
        end
        
        function [p, c] = getSchedulingDist (obj)
            nschedBS =  obj.numUsrSched(:);
            nschedBS(nschedBS == 0) = []; % unscheduled Subframes
            [p, c] = histcounts(nschedBS, 1:obj.nStreams+1);
            p = p/sum(p);            % normalize
            c(end) = [];             % drop the last one
        end
    end
    
    methods (Access='private')
        
        % Scheduling functions are defined with private access. The
        % parameter for scheduling is 
        %               w = SpecEff*QueueLength/ThroughputHistory
        % In a sense this creates a balance between the channel quality,
        % the current data demand and the historical received service. A
        % simple modification to proportional fair with queue length
        % included.
        
        function schedTDMA (obj, ueList, sfnum, bsIndex)
            
            [se, sinr] = obj.DLSinrObj.computeSingleLinkCapacity(ueList);
            queueLen = obj.trafficObj.getQueue(ueList);
            r = obj.thrHistory(ueList)/((sfnum-1)*obj.ttims*1e-3);
            
            wt = (se').*(queueLen > 0)./r;       % Weight for scheduling
            
            [~,ind] = sort(wt,'descend');
            
            cap = obj.eta*obj.bwMHz*(1e6)*se(ind(1)); % in bits

            txSize = min(cap*obj.ttims*1e-3 , queueLen(ind(1)));
            
            obj.trafficObj.dequeue(txSize, ueList(ind(1)));
            
            obj.thrHistory(ueList(ind(1))) = obj.thrHistory(ueList(ind(1))) + txSize;
            obj.totDataTx(ueList(ind(1))) = obj.totDataTx(ueList(ind(1))) + txSize;
            
            obj.fullBufferRate(ueList(ind(1))) = obj.fullBufferRate(ueList(ind(1))) + cap*obj.ttims*1e-3;
            
            obj.numUsrSched(bsIndex, sfnum) = 1;
            
            obj.sinrRecord(ueList) = obj.sinrRecord(ueList) + sinr';
            obj.schedCount(ueList) = obj.schedCount(ueList) + 1;
            
        end
        
        function schedFDMA (obj, ueList, sfnum, bsIndex)
            
            [se, sinr] = obj.DLSinrObj.computeSingleLinkCapacity(ueList);
            queueLen = obj.trafficObj.getQueue(ueList);
            obj.sinrRecord(ueList) = obj.sinrRecord(ueList) + sinr;
            obj.schedCount(ueList) = obj.schedCount(ueList) + 1;
            
            ueList(queueLen == 0) = [];
            se(queueLen == 0) = [];
            queueLen(queueLen == 0) = [];

            r = obj.thrHistory(ueList)/((sfnum-1)*obj.ttims*1e-3);
            
            wt = (se').*(queueLen > 0)./r;       % Weight for scheduling
            
            wt = wt/sum(wt);
            
            bwAlloc = wt*obj.bwMHz;
            
            cap = obj.eta*bwAlloc*(1e6).*(se');
            txSize = min(cap*obj.ttims*1e-3 , queueLen);
            
            obj.fullBufferRate(ueList) = obj.fullBufferRate(ueList) + cap*obj.ttims*1e-3;
            
            obj.trafficObj.dequeue(txSize, ueList);
            
            obj.thrHistory(ueList) = obj.thrHistory(ueList) + txSize;
            
            obj.totDataTx(ueList) = obj.totDataTx(ueList) + txSize;
            
            obj.numUsrSched(bsIndex, sfnum) = length(ueList);
            
            
            
        end
        
        function schedSDMA (obj, ueList, sfnum, bsIndex)
            
            [se, snr] = obj.DLSinrObj.computeSingleLinkCapacity(ueList);
            queueLen = obj.trafficObj.getQueue(ueList);
            
            capPenalty = 0.5;
            
            % remove UEs with empty queues
            se = se'; % Make dimension consistant with the row-major-rest-of-the-code
            se(queueLen == 0)       = [];
            ueList(queueLen == 0)   = [];
            snr(queueLen == 0)      = [];
            queueLen(queueLen == 0) = [];

            r = obj.thrHistory(ueList)/((sfnum-1)*obj.ttims*1e-3);
            wt = se.*(queueLen > 0)./r;       % Weight for scheduling
            [~,ind] = sort(wt,'descend');
            
            if (obj.sdmaICIReject == 1)

                ueSel = ueList(ind(1)); % primary UE scheduled based on best weight
                gam0 = snr(ind(1));
                
                %gammN = 0.2;            % can tuning this lead to better performance?
                pen = 1 - capPenalty;
                gammN = 1/(2^(pen*se(ind(1))) - 1) - 1/gam0; % setting a threshold

                aggresor = ueList(ind(2:end));
                gammaj = obj.DLSinrObj.getIntracellInterference ( ueSel, aggresor);
                
                %updateWt = wt(ind(2:end))./gammaj;
                [~,indx] = sort(gammaj, 'ascend');
                aggresor = aggresor(indx); % sorted in order
            
                gammaj = 10.^(0.1*gammaj); % convert to linear

                % get the point where sum(gamma_j) /gamma_k < Th
                ici = cumsum(gammaj);
                n = ici/gam0;
                ueIndx = sum( n < gammN);

                if (obj.nStreams-1) < ueIndx
                    ueIndx = obj.nStreams - 1;
                end
                
                if (ueIndx > 0)
                    ueSel = [ueSel; aggresor(1:ueIndx)];   
                end
                
            elseif (obj.sdmaICIReject == 2)
                % chooses all the available spatial streams in ascending
                % order of ICI induced. No stopping criteria based on SINR
                % penalty.
                
                ueSel = ueList(ind(1));
                
                aggresor = ueList(ind(2:end));
                gammaj = obj.DLSinrObj.getIntracellInterference ( ueSel, aggresor);
                
                [~,indx] = sort(gammaj, 'ascend');
                aggresor = aggresor(indx);
                
                if (length(indx) < obj.nStreams-1)
                    ueIndx = length(indx);
                else
                    ueIndx = obj.nStreams-1;
                end
                
                if (ueIndx > 0)
                    ueSel = [ueSel; aggresor(1:ueIndx)];   
                end
                
                % does not give better performance for K > 2
                
            else
                %ueIndx = obj.nStreams - 1;
                % take users according to the rate-fairness crteria
                if (length(ueList) >= obj.nStreams)
                    ueSel = ueList(ind(1:obj.nStreams));
                else
                    ueSel = ueList;
                end
            end

            opt.ueSched = ueSel;
            opt.powSplit = length(ueSel);
            
            [specEff, EffSinr] = obj.DLSinrObj.DlSinrCalcBSi(opt);
            
            cap = obj.eta*obj.bwMHz*(1e6).*(specEff'); % full bandwidth is allocated
            queueLen = obj.trafficObj.getQueue(ueSel);
            txSize = min(cap*obj.ttims*1e-3 , queueLen);
            
            obj.fullBufferRate(ueSel) = obj.fullBufferRate(ueSel) + cap*obj.ttims*1e-3;
            
            obj.trafficObj.dequeue(txSize, ueSel);
            obj.thrHistory(ueSel) = obj.thrHistory(ueSel) + txSize;
            obj.totDataTx(ueSel) = obj.totDataTx(ueSel) + txSize;
            obj.numUsrSched(bsIndex, sfnum) = length(ueSel);
            
            obj.sinrRecord(ueSel) = obj.sinrRecord(ueSel) + EffSinr'; % in linear scale
            obj.schedCount(ueSel) = obj.schedCount(ueSel) + 1;
            
        end
        
    end
    
end

