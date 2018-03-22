classdef MimoPlMod < PLMod
    % MimoPlMod: MIMO path loss model
    
    % Public properties
    properties 
        
        % Per cluster path loss factors
        pl0 = 75.85;    % path loss in dB at d=1 meter
        plexp = 3.73;   % Per cluster path loss
        shadStd = 8;    % Log-normal shadowing exponent
        dmin = 1;       % Min distance for path loss computation in meters
        bwH = 7;        % mean horizontal beamwidth in degrees
        
        % number of antennas in horizontal and vertical for both ends 
        nantH = [8 8];   
        nantV = [8 8];
        
        % Other parameters
        ncluster = 3;       % number of clusters                
                                
    end
    
    % Private properties
    properties %(Access=protected)
                
        % To speed computation, the class pre-generates a number of
        % sample covariances
        ncov = 1000;    % number of pre-stored covariances
        cov1;           % per cluster covariance
        cov;            % total covariance  combined over the gains
        sqrtCov;        % sqrt of the covariance combined over the gains
        omniGain;       % omni-directional gain in dB
        bfGainH;        % horizontal BF gain
        bfVecH          % optimal horizontal BF vectors 
        Icov;           % [n1 x n2] matrix of indices for the links
        
        % Cell arrays for each element
        loc;    % (x,y) location 
        
        % Path loss 
        dist;          % Distance in meters between elements
        pathLoss;      % Omn-directional path loss 
        angle;         % LOS angle between elements
        
        % Wrap for computing distances
        % Distances will be computed based on offsets of each column of the
        % matrix wrap
        wrap = [0  0];
    end
    
    methods
        
        % Constructor
        function obj = MimoPlMod()
            obj = obj@PLMod();
        end
        
        % Generates second-order candidate matrices.  
        %
        % To reduce the computation, the method generates a number of sample 
        % covariance matrices which can then be randomly selected and 
        % applied to each link.  We assume a second-order Kronecker MIMO 
        % model of the form:
        %
        %   H = sqrtm(P1)*Q*sqrtm(P2) / Ptot
        %
        % where P1 and P2 are the covariances,
        %   P1 = E[ H*H' ]  P2 = E[ H'*H ]
        %   Ptot = trace(P1) = trace(P2) 
        % and Q is an iid N(0,1) matrix.          
        function genCovMatrices(obj)
                        
            % Get horizontal and vertical angles from the first element
            % in the channel interface vector
            nt  = obj.getNTypes();            
            chanIf1 = obj.getChanIfByInd(1);
            chanIf2 = obj.getChanIfByInd(nt);            
            obj.nantH(1) = chanIf1{1}.nantH;
            obj.nantH(2) = chanIf2{1}.nantH;   
                                   
            % Generate covariance matrices for one cluster            
            ntheta = 100;    % number of points in the theta to integrate over
                        
            % Random beamwidths and central angles
            theta0 = 2*pi*rand(obj.ncov,1);
            bw = -log(rand(obj.ncov,1));
            dtheta = linspace(-1,1,ntheta)';

            % Loop over elemTypes  
            obj.cov1 = cell(2,1);
            for it = 1:nt
                if (obj.nantH(2)==obj.nantH(1)) && (it == 2)
                    obj.cov1{2} = obj.cov1{1};
                    break;
                end
                nant = obj.nantH(it); 
                Q = zeros(nant,nant,obj.ncov);
                for icov = 1:obj.ncov
                    theta = theta0(icov) + bw(icov)*dtheta;
                    h = exp(1i*theta*(1:nant));
                    Q(:,:,icov) = (h'*h)/ntheta;
                end
                obj.cov1{it} = Q;                                
            end
            
           % Generate a random combination
           nant1 = obj.nantH(1);
           nant2 = obj.nantH(nt);
           Q1 = obj.cov1{1};
           Q2 = obj.cov1{nt};
           Qsum1 = zeros(nant1,nant1,obj.ncov);
           Qsum2 = zeros(nant2,nant2,obj.ncov);
           bfVec1 = zeros(nant1,obj.ncov);
           bfVec2 = zeros(nant2,obj.ncov);
           Qsqrt1 = zeros(nant1,nant1,obj.ncov);
           Qsqrt2 = zeros(nant2,nant2,obj.ncov);
           obj.bfGainH = zeros(obj.ncov,1);
           obj.omniGain = zeros(obj.ncov,1);
           for icov = 1:obj.ncov
                % Gemerate shadowing per cluster
                shad = 10.^(0.1*randn(obj.ncluster,1)*obj.shadStd);                
                
                % Select random subset of the per cluster covariances
                I1 = randperm(obj.ncov);
                I2 = randperm(obj.ncov);

                % Add components with random shadowing                
                Qsum1i = zeros(nant1,nant1);
                Qsum2i = zeros(nant2,nant2);
                for ic = 1:obj.ncluster
                    Qsum1i = Qsum1i + Q1(:,:,I1(ic))*shad(ic)*...
                        trace(Q2(:,:,I2(ic)));
                    Qsum2i = Qsum2i + Q2(:,:,I2(ic))*shad(ic)*...
                        trace(Q1(:,:,I1(ic)));
                end
                Qsum1(:,:,icov) = Qsum1i;
                Qsum2(:,:,icov) = Qsum2i;
                Qsqrt1(:,:,icov) = sqrtm(Qsum1i);
                Qsqrt2(:,:,icov) = sqrtm(Qsum2i);
                obj.omniGain(icov) = 10*log10( trace(Qsum1i)/(nant1*nant2));
                
                % Get max BF gain from max singylar vectors
                [U1,S1] = svd(Qsum1i);
                [U2,S2] = svd(Qsum2i);
                bfVec1(:,icov) = U1(:,1);
                bfVec2(:,icov) = U2(:,1);
                obj.bfGainH(icov) = 10*log10(S1(1)*S2(2));
                
           end
           
           % Save results
           obj.cov      = {Qsum1, Qsum2};
           obj.sqrtCov  = {Qsqrt1, Qsqrt2};        
           obj.bfVecH   = {bfVec1, bfVec2};
                       
        end
            
                 
        % Initialization function        
        function init(obj)         
                        
            % Determine number of element types
            nTypes = obj.getNTypes();
            if (nTypes ~= 1) && (nTypes ~= 2)
                error('There must be 1 or 2 element types in the path loss model');
            end         
            
            % Generate the covariance matrices
            obj.genCovMatrices();
            
            % Get the locations of each of the element types
            obj.getLoc();

            % Get the distances between elements
            loc1 = obj.loc{1};
            loc2 = obj.loc{nTypes};
            [obj.dist,obj.angle] = obj.computeDist(loc1,loc2);
                          
            % Compute distance-based path loss
            d = max(obj.dist, obj.dmin);
            obj.pathLoss = obj.pl0 + 10*obj.plexp*log10(d);
            
            % Add antenna gain 
            [n2,n1] = size(obj.pathLoss);     
            chanIf1 = obj.getChanIfByInd(1);
            chanIf2 = obj.getChanIfByInd(nTypes);
            for i1 = 1:n1
                g1 = chanIf1{i1}.antGain(obj.angle(:,i1));
                obj.pathLoss(:,i1) = obj.pathLoss(:,i1) - g1;
            end
            for i2 = 1:n2
                g2 = chanIf2{i2}.antGain(obj.angle(i2,:)');
                obj.pathLoss(i2,:) = obj.pathLoss(i2,:) - g2';
            end
            
            % Generate random indices for each of the links
            obj.Icov = randi(obj.ncov, [n2 n1]);
            
            % Add omni-directional beamforming gain which includes the
            % shadowing
            obj.pathLoss = obj.pathLoss - obj.omniGain(obj.Icov);
        end
                       
        % Compute distance between two vectors of locations
        % loc1 and loc2 should each have 2 columns for the x and y
        % coordinates
        function [dist,angle] = computeDist(obj,loc1, loc2)
            n1 = size(loc1,1);
            n2 = size(loc2,1);
            nwrap = size(obj.wrap,1);
            for iwrap =1:nwrap
                
                % Compute distance at current wrap offset
                offset = obj.wrap(iwrap,:);
                loc1w = loc1 + repmat(offset, n1,1);
                distSq = zeros(n2,n1);
                ndim = size(loc1,2);
                dx = zeros(n2,n1,ndim);
                for idim = 1:ndim
                    dx(:,:,idim) = loc2(:,idim)*ones(1,n1) - ones(n2,1)*loc1w(:,idim)';
                    distSq = distSq + dx(:,:,idim).^2;
                end
                anglei = atan2(dx(:,:,1),dx(:,:,2));
                
                % Check mimima
                if (iwrap == 1)
                    distSqWrap = distSq;
                    angle = anglei;
                else
                    I = (distSq <= distSqWrap);
                    angle = I.*anglei + (1-I).*angle;
                    distSqWrap = min(distSq, distSqWrap);
                    
                end
            end
            dist = sqrt(distSqWrap);            
        end        
        
        % Set wrap for a square area
        function setSqWrap(obj, xmax, ymax)
            % order for wrapping : only works for 2D
                % +
                % 
                % -------------
                % | 1 | 2 | 3 |
                % -------------
                % | 4 | 5 | 6 |
                % -------------
                % | 7 | 8 | 9 |
                % -------------    +
            if (xmax~=0)&&(ymax~=0)
                obj.wrap = [-xmax ymax; 0 ymax; xmax ymax; -xmax 0;...
                            0 0; xmax 0; -xmax -ymax; 0 -ymax; xmax -ymax];
            end
        end
        
        % Get the wrap matrix
        function wrap = getWrap(obj)
            wrap = obj.wrap;
        end
        
        % Get the path loss matrix
        function pathLoss = getPathLoss(obj, name1, name2, ind1, ind2)
            [exist,rev] = obj.getOrder(name1,name2);
            if ~exist                
                error('No path loss between %s and %s', name1, name2);
            end
            if (rev)
                pathLoss = obj.pathLoss';
            else
                pathLoss = obj.pathLoss;               
            end
            if (nargin >= 4)
                if (~isempty(ind1))
                    pathLoss = pathLoss(:,ind1);
                end
            end
            if (nargin >= 5)
                if (~isempty(ind2))
                    pathLoss = pathLoss(ind2,:);
                end
            end
            
        end        
        
        % Gets the locations of the mobiles
        function loc = getLoc(obj)
            nTypes = obj.getNTypes();
            obj.loc = cell(nTypes,1);            
            for itype= 1:nTypes
                chanIfVec = obj.getChanIfByInd(itype);
                nelem = length(chanIfVec);
                loci = zeros(nelem,2);
                for ielem = 1:nelem
                    loci(ielem,:) = ...
                        chanIfVec{ielem}.getDev().getNode().getMobMod().getLoc()';                    
                end
                obj.loc{itype} = loci;
            end
            loc = obj.loc;
        end
        
        function dist = getDist(obj)
            dist = obj.dist;
        end
        
    end
    
end

