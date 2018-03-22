classdef UnifLayout < hgsetget
    % UnifLayout:  Layouts 
    properties (Access = private)        
        % Dimensions of the rectangle
        xmin, xmax, ymin, ymax;                
    end
    
    methods
        
        % Constructor
        function obj = UnifLayout(posLim)
            obj.xmin = 0;
            obj.xmax = posLim(1);
            obj.ymin = 0;
            obj.ymax = posLim(2);
        end
        
        % Randomly drop locations in the locations
        function loc = drop(obj, nNodes)
            loc = rand(nNodes,2);
            loc(:,1) = loc(:,1)*(obj.xmax-obj.xmin)+obj.xmin;
            loc(:,2) = loc(:,2)*(obj.ymax-obj.ymin)+obj.ymin;            
        end
    end
    
    
end