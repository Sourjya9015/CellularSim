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
        psign; % signal power dB
        pnoise; % noise power dB
        pintrfr; % interference power dB
                
        % Number of UEs per base station
        nueBS;
        
        % Bandwidth in MHz assigned to each UE
        bwMHz;
                
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
            bwMHzTot = opt.bwMHzTot;    % bandwidth in MHz
            Icell = opt.Icell;            % Icell(j) = RX index for TX i
            [nue,nbs] = size(opt.pathLoss); % Dimensions
            
            multiacc = 'tdma'; % default
            kstreams = 1; % default
            if nargin == 2
                multiacc = multiacsOpt.multiacs;
                kstreams = multiacsOpt.K;
                if (kstreams <= 0)
                   error('Invalid value in  DLSinrCalc: multiacsOpt.K <= 0');
                end
            end

            txpow = txpow(:)';

            
            % SD: Tx power need to be divided into K streams
            % for SDMA
            if (strcmpi(multiacc,'sdma'))
                %disp('sdma transmission!');
                txpow = txpow - 10*log10(kstreams);
            end
            
            p = repmat(txpow,nue,1) - pathLoss;
            p = 10.^(0.1*p);
            
            sigpow = zeros(nue,1);
            for iue = 1:nue
                sigpow(iue) = p(iue,Icell(iue));
            end
            
            obj.nueBS = zeros(nbs,1);
            for iue = 1 : nue
                obj.nueBS(Icell(iue)) = obj.nueBS(Icell(iue)) + 1;
            end
            
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
            
            % SINR & INR calculation
            noisepow = 10.^(0.1*noisepow);
            totPow = sum(p,2) + noisepow;
            sigPowdB = 10*log10(sigpow);
            
            iplusn = 10*log10(totPow - sigpow);  % interference and noise
            obj.sinr = sigPowdB - iplusn;       % SINR
            obj.inr = 10*log10(totPow - sigpow - noisepow) - 10*log10(noisepow); % interference to noise ratio
            
            %nb:
            %powers in linear:
            obj.psign = sigpow;
            obj.pnoise = noisepow;
            obj.pintrfr = totPow - sigpow - noisepow;
            
            
            % rate calculation, Mbps/cell (@BW = obj.bwMHz)
            obj.bwMHz = zeros(nue,1);
            for iue = 1 : nue
                obj.bwMHz(iue) = kstreams*bwMHzTot/obj.nueBS(Icell(iue)); % cell capacity
            end            
            beta = 10^(-0.1*obj.sirLossdB);
            se = log2(1 + beta*(10.^(0.1*(obj.sinr)))); % spectral efficiency
            obj.specEff = se;
            obj.rate = obj.bwMHz.*min(se, obj.maxSe);            
        end
        
    end
    
end

