classdef HexLayout
    properties
        isd;        % inter-site distance
        nx, ny;     % num of cols and rows of hexagon
    end
    
    methods
        
        % Constructor.
        % HexLayout(isd,nx,ny) creates a hex layout with nx x ny cell sites
        % HexLayout(isd,posLim) computes nx and ny to roughly match the
        % position limts.       
        function obj = HexLayout(isd, arg1, arg2)
            obj.isd = isd;
            
            % If only one argument is specified, then it will automatically
            % nx and ny from the posLim
            if (nargin <= 2)
                posLim = arg1;
                dy = (isd/2)*(1 + 1/sqrt(3));
                obj.nx = round(posLim(1)/obj.isd);
                obj.ny = round(posLim(2)/dy);
            else
                obj.nx = arg1;
                obj.ny = arg2;
            end
            
        end
        
        % Generate locations for centers of the hexagons
        function loc = drop(obj,nnodes)
            
            % Check if the number of nodes is correct
            if (nnodes ~= obj.nx*obj.ny)
                error('Incorrect number of nodes');
            end
            
            % Create node positions
            dy = 1 + 1/sqrt(3);
            x = repmat((0:2:(2*obj.nx-2)),obj.ny,1);
            x((1:2:obj.ny),:) = x((1:2:obj.ny),:) + 1;
            y = repmat((0:obj.ny-1)'*dy,1,obj.nx);
            x = (obj.isd/2)*x(:);
            y = (obj.isd/2)*y(:);
            loc = [x y];
        end
        
        % Generate position limits
        function posLim = getPosLim(obj)
            dy = 1 + 1/sqrt(3);
            posLim = [2*obj.nx dy*obj.ny]*obj.isd/2;
        end
        
        % Get number of elements
        function [nx,ny] = getDim(obj)
            nx = obj.nx;
            ny = obj.ny;
        end
    
    end
end

