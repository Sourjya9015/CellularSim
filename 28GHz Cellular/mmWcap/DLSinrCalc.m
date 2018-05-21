classdef DLSinrCalc
    %DLSINRCALC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % sinr
        sinr;
        % inr
        inr;
        % rate
        rate;
        % spectral efficiency
        specEff;
        
        %nb
        psign;      % signal power dB
        pnoise;     % noise power dB
        pintrfr;    % interference power dB
        alpha = 0;      % inverse coding gain due to quant.
        rxBFgains;  % Rx BF gains (dB)
        
        psdNoise; % Noise power spectral density.
                
        % Number of UEs per base station
        nueBS;
        
        % Bandwidth in MHz assigned to each UE
        bwMHz;
        bwMHzTot;
        maxStreams = 1; % Only set for SDMA systems
        intraPow;
                
        % Options
        sirLossdB = 3;      % Loss from Shannon capacity
        maxSe = 4.8;        % max spectral efficiency in bps/Hz
    end
    
    methods
        function obj = DLSinrCalc(opt, multiacsOpt)
            % Get parameters
            pathLoss = opt.pathLoss;    % path loss in dB  (nue X nbs)
            txpow = opt.txpow;          % max TX pow in dBm
            noisepow = opt.noisepow;    % noise power in dBm   
            obj.bwMHzTot = opt.bwMHzTot;  % bandwidth in MHz
            Icell = opt.Icell;            % Icell(j) = RX index for TX i
            [nue,nbs] = size(opt.pathLoss); % Dimensions
            
            intraPL = opt.intraOpt;    % Only needed for SDMA systems; 
            % No intra-cell intrefernce when orthogonal Tx is used, i.e
            % (TDMA/OFDMA)
            
            multiacc = 'tdma'; % default
            kstreams = 1; % default
            if nargin == 2
                multiacc = multiacsOpt.multiacs;
                
            end

            txpow = txpow(:)';

            
            % Tx power divided into K streams for SDMA
            if (strcmpi(multiacc,'sdma'))
                %disp('sdma transmission!');
                kstreams = multiacsOpt.K;
                txpow = txpow - 10*log10(kstreams);

                obj.maxStreams = kstreams;
                if (kstreams <= 0)
                   error('Invalid value in  DLSinrCalc: multiacsOpt.K <= 0');
                end
                
            end
            
            p = repmat(txpow,nue,1) - pathLoss;
            p = 10.^(0.1*p);
            
            sigpow = zeros(nue,1);
            for iue = 1:nue
                sigpow(iue) = p(iue,Icell(iue));    % Tx power/pathLoss (lin scale)
            end
            
            obj.nueBS = zeros(nbs,1);
            obj.intraPow = containers.Map();
            intraIntfPow = zeros(nue,1);
            
            for iue = 1 : nue
                obj.nueBS(Icell(iue)) = obj.nueBS(Icell(iue)) + 1;
                
                % for SDMA
                if (strcmpi(multiacc,'sdma'))
                    txPower = txpow(Icell(iue));
                    intraCellIntf = intraPL(num2str(iue));
                    intraCellInfo.ueList = intraCellIntf.ueList;
                    intraCellInfo.intrfpow = (txPower + 10*log10(kstreams)) - intraCellIntf.pathLoss;
                    obj.intraPow(num2str(iue)) = intraCellInfo;
                    
                    intraCellin = 10.^(0.1*(txPower - intraCellIntf.pathLoss));
                    
                    intraCellin = sort(intraCellin);

                    intraIntfPow(iue) = sum(intraCellin(1:kstreams-1)); % in an average sense
                    % One of the K streams are Tx towards the UE.
                end
            end
            
            %intraIntfPow = 10*log10(intraIntfPow);
            
            k = 0;
            disabledTx = zeros(nbs,1);
            for ibs = 1 : nbs
                if obj.nueBS(ibs)==0
                    k = k+1;
                    disabledTx(k) = ibs;
                end
            end
            disabledTx = disabledTx(1:k);
            
            p(:,disabledTx) = 0;
            
            % rate calculation, Mbps/cell (@BW = obj.bwMHz)
            obj.bwMHz = zeros(nue,1);
            for iue = 1 : nue
                obj.bwMHz(iue) = kstreams*obj.bwMHzTot/obj.nueBS(Icell(iue)); % cell capacity
            end
            
            obj.psdNoise = noisepow - 10*log10(obj.bwMHzTot); % Make sure this has a +60 for the MHz conversoin in the main script
            
            if(length(noisepow) == 1)
                noisepow = repmat(noisepow,nue,1);
            end
            
            % SINR & INR calculation
            noisepow = 10.^(0.1*noisepow);
            totPow = sum(p,2) + noisepow + intraIntfPow;

            sigPowdB = 10*log10(sigpow);
            %iplusn = 10*log10(totPow - sigpow);  % interference and noise
            
            %obj.sinr = sigPowdB - iplusn;       % SINR
            
            %powers in linear:
            obj.psign = sigpow;
            obj.pnoise = noisepow;
            obj.pintrfr = totPow - sigpow - noisepow - intraIntfPow; % The inter Cell Intf only

            iplusn = 10*log10(noisepow + obj.pintrfr + intraIntfPow);  % interference and noise
           
            %nb : calculate the quantized SINR
            obj.alpha = opt.alpha;
            obj.rxBFgains = opt.rxBFgains;
            gamma = 10.^(0.1*(sigPowdB - iplusn));              % pre-quantization SINR (lin)
            
            % compute the post quantization SINR
            if (obj.alpha > 0)
                n0 = noisepow + obj.pintrfr;
                g0 = sigPowdB - n0;
                gj = intraIntfPow - n0;
                
                g0_preBF = g0 ./ (10.^(0.1*obj.rxBFgains));   
                gj_preBF = gj ./ (10.^(0.1*obj.rxBFgains));  
                
                gamma = ( (1 - obj.alpha)*g0 )./(1 + (1-obj.alpha)*gj+ obj.alpha*(g0_preBF + gj_preBF));       % SINR after quantization
            
            end
            
            obj.sinr = 10*log10(gamma);
            obj.inr = 10*log10(totPow - sigpow - noisepow) - 10*log10(noisepow); % interference to noise ratio

            %nb:           
            beta = 10^(-0.1*obj.sirLossdB);
            se = log2(1 + beta*(10.^(0.1*(obj.sinr)))); % spectral efficiency
            obj.specEff = se;
            obj.rate = obj.bwMHz.*min(se, obj.maxSe);            
        end
        
        % Computes the SINR and SE for each BS. The BW (opt.bw) should be
        % set using the particular scheduling algorithm.
        % return the SINR and the Spectral efficiencies for data Tx.
        
        % nb: This fn. works properly only if DlSinrCalc has been called
        % before. It is a weird dependency!
        % sd: No, DlSinrCalc is not a function, it is the constructor of
        % the class. It will be called for any and every object
        % instantiation. The dependency is natural and useful.
        
        function specEff = DlSinrCalcBSi (obj, opt) 
            
            ueSched = opt.ueSched;     % List of UEs scheduled in this instant
            %bwSched = opt.bwSched;    % a vector of bandwidth allocated to each user
            splitK = opt.powSplit;
            
            noisePow = obj.psdNoise + 10*log10(obj.bwMHzTot);
            
            % already store in Linear scale
            sigPow = obj.psign(ueSched);
            % do appropriate power division when max streams are not used
            if (obj.maxStreams > splitK)
                sigPow = sigPow*obj.maxStreams/splitK; % in linear scale   
            end
            
            % SDMA case; intra cell interference if the intacell Map is
            % filled in applyBfGain.
            
            intraIntf = zeros(length(ueSched),1);
            if (splitK > 1 && ~isempty(obj.intraPow))
                
                for iue = 1:length(ueSched)
                    
                    othUsr = ueSched;
                    othUsr(othUsr == ueSched(iue)) = [];
                    intraInfo = obj.intraPow(num2str(ueSched(iue)));
                    
                    if isempty(intraInfo)
                        error('Programming error: check how the map is accessed');
                    end
                    
                    [~,indx,~] = intersect(intraInfo.ueList, othUsr);
                    intraIntf(iue) = sum(10.^(0.1*(intraInfo.intrfpow(indx))))/splitK;
                    
                end
                
            end

            intf = obj.pintrfr(ueSched);
            
            noisePow = 10.^(0.1*noisePow);
            % Intra Cell Intf in case of SDMA

            %nb:
            gamma_0 = sigPow./(intf + noisePow );
            gamma_ici = intraIntf./(intf + noisePow);
            gamma_0sansRxBF = gamma_0 ./ (10.^(0.1*obj.rxBFgains(ueSched)));
            gamma_icisansRxBF = gamma_ici ./ (10.^(0.1*obj.rxBFgains(ueSched)));
            % Sinr after quantization with intra cell and rx beamforming
            Sinr = ((1 - obj.alpha)*gamma_0)./(1 +  (1-obj.alpha)*gamma_ici + obj.alpha*(gamma_0sansRxBF + gamma_icisansRxBF)); 
            
            beta = 10^(-0.1*obj.sirLossdB);
            specEff = log2(1 + beta*Sinr);
            specEff = min(specEff, obj.maxSe);
        end
        
        % Used for scheduling to get a metric for the channel quality.
        function [se, gamma] = computeSingleLinkCapacity (obj, ueList)
            
            noisePow = obj.psdNoise + 10*log10(obj.bwMHzTot);   %  +60 is included in the main script
            noisePow = 10.^(0.1*noisePow); % linear scale
            
            sigPow = obj.psign(ueList)*obj.maxStreams;
            
            intf = obj.pintrfr(ueList);
            
            % take care of quantization penalty
            gamma = sigPow./(intf + noisePow); % returns the SINR without quantization
            gamma_preRx = gamma ./ (10.^(0.1*obj.rxBFgains(ueList)));
            Sinr = ((1 - obj.alpha)*gamma)./(1 + obj.alpha*gamma_preRx);
            
            beta = 10^(-0.1*obj.sirLossdB);
            se = log2(1 + beta*Sinr);
            se = min(se, obj.maxSe);
            % se in bits/sec/Hz, gamma in linear scale
        end
        
        function gammaj = getIntracellInterference (obj, victim, aggressors)
            % returns (P*BF_Gain)/(N0 + Z0) for all UEs in the same cell a the
            % victim UE. The full power is  multiplied and needs to be modified by
            % P/K during scheduling when the total number of simultaneous
            % non-othogonal transmissions are decided. Full bandwidth
            % allocation is assumed for each user.
            
            if (isempty(obj.intraPow))
                error('Programming error: Interference power not calculated');
            end
            intraInfo = obj.intraPow(num2str(victim));
            
            
            % the UE list should have all the aggressors; the reverse case
            % need not always hold.
            [~,indxa, indxb] = intersect (intraInfo.ueList, aggressors);
            
            iciPow = zeros(1,length(aggressors));
            iciPow(indxb) = intraInfo.intrfpow(indxa); % arranged in order
            
            % convert the received Intra-Cell Intf powers to ICI SINRs
            gammaj = 10.^(0.1*iciPow)/(obj.pnoise(victim) + obj.pintrfr(victim));
            gammaj = 10*log10(gammaj);
            % returns in dB
            
        end
        
    end
    
end

