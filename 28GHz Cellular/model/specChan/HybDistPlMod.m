classdef HybDistPlMod < DistPlMod
    % DistPlMod Distance-based path loss model with shadowing.  
    %
    % Assumes a path loss model of the form:
    %
    %   PL(d) = pl0 + plexp*10*log10(d) + shad
    %
    % where pl0 is the average path loss at 1 meter and stad is 
    % lognormal shadowing
    % 
    % (pl0, plexp) is chosen as (pl0LOS, plexpLOS) or (pl0NLOS, plexpNLOS)
    % with probability p(d) and 1 - p(d)
    
    % Public properties
    properties 
        
        plexpLOS = 3.7; % Path loss exponent for Line of Sight
        plexpNLOS = 3.7;% Path loss exponent for NonLine of Sight
        plexpOut = 10; % A very high path loss exponent for outage
        
        % The constant term in the path loss model
        % If pl0 is empty, then it is automatically computed to
        % match the freespace propagation from Friis formula.
        pl0LOS = 10;    % Path loss in dB at d=1 meter for Line of Sight
        pl0NLOS = 10;   % Path loss in dB at d=1 meter for NonLine of Sight
        pl0Out = 300;   % A very high Path loss in dB for outage case
        
        % fuction handler for LOS probability based on distance
        probR = @(R)min(18./R,1).*(1-exp(-R/63))+exp(-R./63);
        % function handler for outage probability based on distance
        poutR = @(R)zeros(size(R));
        
        % a link can be either:
        % outage,   with probability poutR, or
        % LOS,      with probability (1-poutR)*probR, or
        % NLOS      with probability (1-poutR)*(1-probR)
        
        
        shadStdLOS;
        shadStdNLOS;
        
        region;         % matrix of logicals, true if UEi is in LOS of BSj
        regionout;      % matrix of logicals, true if UEi is out of sight of BSj
    end
    
%     % Private properties
%     properties (Access=protected)
%         
%     end
    
    methods
        % Initialization function
        % Set the shadowing
        function init(obj)         
            
            % Compute the constant term, if requested
            if (isempty(obj.pl0NLOS))
                obj.setPl0NLOS();
            end
            if (isempty(obj.pl0LOS))
                obj.setPl0LOS();
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
             
            % Compute regions: outage, LOS or NLOS
            obj.calcRegion();
            
            % Shadowing 
            obj.genShad();
            
            % Compute path loss
            d = max(obj.dist, obj.dmin);
            pathLossLOS = obj.pl0LOS + 10*obj.plexpLOS*log10(d);
            pathLossNLOS = obj.pl0NLOS + 10*obj.plexpNLOS*log10(d);
            pathLossOut = obj.pl0Out + 10*obj.plexpOut*log10(d);
            obj.pathLoss = (1-obj.regionout).*obj.region.*pathLossLOS +...
                (1-obj.regionout).*(1-obj.region).*pathLossNLOS +...
                obj.regionout.*pathLossOut + obj.shad;
            
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
        
       
        % Computes the constant term in the path loss to match 
        % Friis' freespace path loss:
        %   PL(dref) = -20*log10(lam/4/pi/dref) = pl0NLOS + 10*plexpNLOS*log10(dref)
        function setPl0NLOS(obj)
            c   = 3e8;                  % velocity of light in m/s
            lam = c/(1e6*obj.fcMHz);    % wavelength
            obj.pl0NLOS = -20*log10(lam/4/pi/obj.dref) - ...
                10*obj.plexpNLOS*log10(obj.dref);
        end
        
        function setPl0LOS(obj)
            c   = 3e8;                  % velocity of light in m/s
            lam = c/(1e6*obj.fcMHz);    % wavelength
            obj.pl0LOS = -20*log10(lam/4/pi/obj.dref) - ...
                10*obj.plexpLOS*log10(obj.dref);
        end

        
        function calcRegion(obj)
            p = feval(obj.probR,obj.dist);
            r = rand(size(p));
            obj.region = r<p;
            
            p = feval(obj.poutR,obj.dist);
            r = rand(size(p));
            obj.regionout = r<p;
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
            shadLOS = randn(n2,n1)*obj.shadStdLOS;
            shadNLOS = randn(n2,n1)*obj.shadStdNLOS;
            obj.shad = obj.region.*shadLOS + (1-obj.region).*shadNLOS;
        end
        
    end
    
end

