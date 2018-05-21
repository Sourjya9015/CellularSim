classdef trafficGen < hgsetget
    % trafficGen:  Generates traffic patterns
    %  Simulates a poisson arrival process, with exponntial packet sizes
    
    properties (Access=private)
        % Always generates Poisson processes
        lambda;   % arrival rates, pkt/sec
        pktsize;  % mean packet size, in bits
        frac;     % the fraction of users belonging to each flow
        
        nflows;
        nqueue;   % number of queues
        queues;   % data queues
        arrRates; % per user arrival rate
        
        totDat;   % total data arriving in the queue during the simulation epoch 
        totTime;
        cumQueueLen; % cumulative queue length, used to calculate the average
        arrivals; % arrival stat,arrivals.T = arrival time(1xN) and arrivals.type = maps to data dist.
        
        tti;      % transmission time intervals
        
    end
    
    methods
        
        function obj = trafficGen(trafficOpt)
            
            obj.lambda = trafficOpt.lambda;
            obj.pktsize = trafficOpt.size*8; % convert bytes to bits
            obj.frac = trafficOpt.frac;
            
            if (length(obj.lambda) ~= length(obj.pktsize) ...
                    && length(obj.lambda) ~= length(obj.frac)...
                    && ~isempty(obj.lambda)~=0 && ~isempty(obj.pktsize)~=0)
               error('Invalid parameters, arrival rates and packet sizes should have same length');
            end

            if (sum(obj.frac)~=1)
               error('Parameter frac in  trafficGen is incorrectly set');
            end
            obj.nflows = length(obj.lambda);
            
            obj.nqueue = trafficOpt.nqueue;
            obj.queues = zeros(1,obj.nqueue);
            obj.cumQueueLen = zeros(1,obj.nqueue);
            obj.arrRates = zeros(1,obj.nqueue);
            
            obj.tti = trafficOpt.tti; % the Transmission time intervals
            obj.totTime = ceil(trafficOpt.totTime); % total time in seconds
            
            numSF = obj.totTime/obj.tti;
            % generate arrival times
            
            inds = randperm(obj.nqueue);
            %disp(inds);
            k = 1;
            for icnt = 1:length(obj.lambda)
                
                nflow = round(obj.frac(icnt)*obj.nqueue);
                
                
                for jcnt = k:(k+nflow-1)
                    arv = false(1,numSF); % saves size, use logicals
                    rnd = rand(size(arv));
                    arv(rnd < obj.lambda(icnt)*obj.tti) = true;    % sets the SFs where arrivals occur
                    obj.arrivals(inds(jcnt)).T    = arv;  % 0-1 list of SFs with arrivals/no arrival
                    obj.arrivals(inds(jcnt)).type = icnt; % Maps to the size of packet
                end
                
                % convert to average b/s from pkt/s: lambda (pkts/s) * avgPktSize (bits/pkts)
                obj.arrRates(inds(k:(k+nflow-1))) = obj.lambda(icnt)*obj.pktsize(icnt);
                k = k+nflow;
            end
            
            obj.totDat = zeros(1,obj.nqueue);
        end
        
        function tti = getTTI(obj)
           tti =  obj.tti; 
        end
        
        function dequeue(obj,nbits,ind)
            
            if (max(ind) > obj.nqueue)
               error('Index out of bound for traffic queue'); 
            end
            
            obj.queues(ind) = max(obj.queues(ind) - nbits, 0);    
        end
        
        
        
        function [size] = getQueue(obj, ind)
            size = obj.queues(ind);
        end
        
        function genTraffic(obj, sfnum)
            
            % Update the avarage queue length
            % compute queue length here. Packets arriving now are
            % served in the next SF at the earliest.
            obj.cumQueueLen = obj.cumQueueLen+ obj.queues;
            
            % Check, can this be vectorized?
            for usr=1:obj.nqueue

                if(obj.arrivals(usr).T(sfnum) == true)
                    muPkt = obj.pktsize(obj.arrivals(usr).type);
                    
                    pktSz = exprnd(muPkt);  % exponentially distributed packet sizes
                    obj.totDat(usr) = obj.totDat(usr) + pktSz;
                    
                    obj.queues(usr) = obj.queues(usr) + pktSz;
                end 
            end
            
            
        end

        
        function totQueue = getTotDataQueues(obj)
            totQueue = obj.totDat;
        end
        
        function [Eq, T] = calcAvgQueueLength (obj, numSFtot)
            Eq = obj.cumQueueLen/numSFtot;
            
            T = Eq./obj.arrRates;
        end
    end
    
end

