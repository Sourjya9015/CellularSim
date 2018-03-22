classdef DistPlMod < PLMod
    % DistPlMod Distance-based path loss model with shadowing.  
    %
    % Assumes a path loss model of the form:
    %
    %   PL(d) = pl0 + plexp*10*log10(d) + shad
    %
    % where pl0 is the average path loss at 1 meter and stad is 
    % lognormal shadowing
    
    % Public properties
    properties 
        plexp = 3.7;    % Path loss exponent
        shadStd = 8;    % Log-normal shadowing exponent
        dmin = 1;       % Min distance for path loss computation in meters
        
        % The constant term in the path loss model
        % If pl0 is empty, then it is automatically computed to
        % match the freespace propagation from Friis formula.
        pl0 = 10;       % Path loss in dB at d=1 meter
        dref = 3;       % Reference distance in meters
                        % The intercept is calculated to meet the 
        fcMHz = 2100;   % Carrier freq in MHz
    end
    
    % Private properties
    properties (Access=protected)
        
        % Cell array of names and locations for each element type
        loc;
        
        % Matrices stored by 
        dist;   % Distance in meters between elements
        shad;   % Shadowing in dB between elements
        pathLoss;   % Path loss in dB 
        angle;  % LOS angle between elements
        
        % Wrap for computing distances
        % Distances will be computed based on offsets of each column of the
        % matrix wrap
        wrap = [0  0];
    end
    
    methods
        % Initialization function
        % Set the shadowing
        function init(obj)         
            
            % Compute the constant term, if requested
            if (isempty(obj.pl0))
                obj.setPl0();
            end
            
            % Determine number of element types
            nTypes = obj.getNTypes();
            if (nTypes ~= 1) && (nTypes ~= 2)
                error('There must be 1 or 2 element types in the path loss model');
            end                        
            
            % Get the locations of each of the element types
            obj.getLoc();

            % Get the distances between elements
            loc1 = obj.loc{1};
            loc2 = obj.loc{min(2,nTypes)};
            [obj.dist,obj.angle] = obj.computeDist(loc1,loc2);
             
            % Shadowing 
            obj.genShad();
            
            % Compute path loss
            d = max(obj.dist, obj.dmin);
            obj.pathLoss = obj.pl0 + 10*obj.plexp*log10(d) + obj.shad;
            
            % Add antenna gain -- will add momentarily
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
        end
        
        % Generates a random shadowing matrix, obj.shad.
        % If nTypes == 2, then obj.shad(i2,i1) is the
        % shadowing between the i1-th element of the 1st interface type and
        % i2-th element of the 2nd interface type.  If nTypes == 1, then
        % the obj.shad(i2,i1) is the shadowing between the i1 and i2-th
        % element of the 1st interface.
        function genShad(obj)
            
            % Get dimensions
            nTypes = obj.getNTypes();            
            n1 = length( obj.getChanIfByInd(1) );
            if (nTypes == 2)
                n2 = length( obj.getChanIfByInd(2) );
            else
                n2 = n1;
            end
            
            % Generate random shadowing
            obj.shad = randn(n2,n1)*obj.shadStd;                        
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
            obj.wrap = [-xmax ymax; 0 ymax; xmax ymax; -xmax 0; 0 0; xmax 0; -xmax -ymax; 0 -ymax; xmax -ymax];
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
        
        % Computes the constant term in the path loss to match 
        % Friis' freespace path loss:
        %   PL(dref) = -20*log10(lam/4/pi/dref) = pl0NLOS + 10*plexpNLOS*log10(dref)
        function setPl0(obj)
            c   = 3e8;                  % velocity of light in m/s
            lam = c/(1e6*obj.fcMHz);    % wavelength
            obj.pl0 = -20*log10(lam/4/pi/obj.dref) - ...
                10*obj.plexp*log10(obj.dref);
        end

        
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

