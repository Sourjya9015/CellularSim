classdef AntPattern < hgsetget
    % AntPattern:  Base class for an antenna pattern    
    properties
        sectInd;    % sector index
    end
    methods 
        % Constructor is based on the sector index
        function obj = AntPattern(sectInd)
            if (nargin < 1)
                sectInd = 1;
            end
            obj.sectInd = sectInd;            
        end
        % Returns gain as a function of the angle in radians
        function gaindB = gain(obj, angle)
            gaindB = 0;
        end
    end
    
end

