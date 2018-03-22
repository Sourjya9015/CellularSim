classdef ULPowConOptFn < hgsetget
    % ULPowConOptFn:  Optimization function for uplink power control
    properties
        % Path loss matrices such that the SINR on link i is
        % p(i)/z(i), z = G*p+g0 where
        % p is the vector of transmit power fractions for transmitter i
        G, g0;
        
        % Number of UEs per base station
        nueBS;
        
        % Bandwidth in MHz assigned to each UE
        bwMHz;
        
        % Total bandwidth
        bwMHzTot;
        
        % Boolean flag indicating if FDMA is enabled
        fdma;
        
        % Options
        sirLossdB = 3;      % Loss from Shannon capacity
        maxSe = 4.8;        % max spectral efficiency in bps/Hz
        
        % Optimization options
        verbose = false; % print optimization progress
        nit = 10;        % num iterations
        powMin = 0.01;   % min power in linear scale relative to max        
    end
    
    methods
        
        % Contructor
        function obj = ULPowConOptFn(opt)
            
            % Get parameters
            pathLoss = opt.pathLoss;    % path loss in dB
            txpow = opt.txpow;          % max TX pow in dBm
            noisepow = opt.noisepow;    % noise power in dBm
            obj.bwMHzTot = opt.bwMHzTot;    % bandwidth in MHz
            Isel = opt.Isel;            % Isel(j) = RX index for TX i
            [nbs,nue] = size(opt.pathLoss); % Dimensions
            obj.fdma = opt.fdma;            % mulitple access mode
            
            
            % Compute gain matrix G and g0
            obj.g0 = zeros(nue,1);
            obj.G = zeros(nue,nue);
            for iue = 1:nue
                ibs = Isel(iue);
                obj.G(iue,:) = txpow'-pathLoss(ibs,:);
                obj.g0(iue) =  noisepow-obj.G(iue,iue);
                obj.G(iue,:) = obj.G(iue,:)-obj.G(iue,iue)- 10*log10(sum(Isel==Isel(iue)));
                obj.G(iue,iue) = 0;
            end
            obj.G = 10.^(0.1*obj.G);
            obj.g0 = 10.^(0.1*obj.g0);
            
            % Count number of UE connected to each BS and zero out the
            % interference from UEs connected to the same BS. Also
            % set the bw for each UE
            obj.nueBS = zeros(nbs,1);
            obj.bwMHz = zeros(nue,1);
            for ibs = 1:nbs
                I = find(Isel == ibs);
                obj.nueBS(ibs) = length(I);
                obj.bwMHz(I) = obj.bwMHzTot/obj.nueBS(ibs);
                for iue = I
                    obj.G(iue,I) = 0;
                end
                
                % If FDMA is used, then the SNR is boosted by the bandwidth
                % fraction
                if (obj.fdma)
                    obj.G(I,:) = obj.G(I,:) / obj.nueBS(ibs);
                    obj.g0(I) = obj.g0(I) / obj.nueBS(ibs);
                end
            end
        end
        
        % Measure SINR for a power vector p
        function sinr = getSinr(obj, p)
            if (nargin < 2)
                nue = length(obj.g0);
                p = ones(nue,1);
            end
            sinr = 10*log10(1./(obj.G*p + obj.g0));
            
        end
        
        % Computes spectral efficiency as a function of the transmit power
        % vector.  Also computes the derivitive of log-sum utility with
        % respect to the power.
        function [se,dutil,sinr,z0] = computeSE(obj, p)
            
            % Compute interference
            z0 = obj.G*p;
            z = z0 + obj.g0;
            
            % Compute SE
            sinr = p./z;
            beta = 10^(-0.1*obj.sirLossdB);
            se = log2(1 + beta*sinr);
            I = (se <= obj.maxSe);
            se = obj.bwMHz.*min(se, obj.maxSe);
            
            % Compute derivitive of dratew = w'*(drate / dp)
            if (nargout >= 2)
                
                % dse = dse / dsinr
                dse = beta/log(2)*obj.bwMHz ./(1+beta*sinr);
                dutil = dse.*I./se;
                dutil = dutil.*(1./z) - obj.G'*(dutil.*(p./(z.^2)));
            end
            
        end
        
        % Compute rate, SINR and INR given a power allocation
        function [rate, sinr, inr] = computeRate(obj,p)
            
            % Get dimensions
            nue = length(obj.g0);
            if (length(p) == 1)
                p = repmat(p,nue,1);
            end
            
            % Compute SINR
            [rate,~,sinr,z0] = computeSE(obj, p);
            
            % Compute INR
            inr = z0./obj.g0;
            
        end
        
        % Optimize power
        function p = optPower(obj)
            
            % Initial power at maximum 
            nue = length(obj.g0);
            p = ones(nue,1);
            
            % Initial rate and gradient
            [se,utilGrad] = obj.computeSE(p);
            util = sum(log(se));
            step = 0.01;
            
            for it = 1:obj.nit
                
                % Compute candidate point
                p1 = p + step*utilGrad;
                p1 = max(obj.powMin,p1);
                p1 = min(p1, 1);
                
                % Get new utility and gradient
                [se1,utilGrad1] = obj.computeSE(p1);
                util1 = sum(log(se1));
                
                % Check if pass
                alpha = 0.5;
                pass = (util1-util > alpha*utilGrad'*(p1-p));
                
                % Update candidate and step
                if pass
                    p = p1;
                    utilGrad = utilGrad1;
                    util = util1;
                    step = step*2;
                else
                    step = step*0.5;
                end
                
                % Print progress
                if (obj.verbose)
                fprintf(1,'it=%d util=%12.4e step=%12.4e\n', it, ...
                    exp(util/nue), step);
                end
            end
        end
        
    end
    
end

