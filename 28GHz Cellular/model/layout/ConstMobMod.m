classdef ConstMobMod < MobMod
    % ConstMobMod:  Constant mobility model
    properties
        loc = [];    % 2x1 or 3x1 vector for position.  
    end
    
    methods
        % Constructor
        function obj = CartPosMod(loc)  
            if (nargin >= 1)
                obj.loc = loc;            
            end
        end
        
        % Accessor
        function [loc] = getLoc(obj)
            loc     = obj.loc;
        end
        function setLoc(obj, loc)
            obj.loc = loc;            
        end
        
    end
    
end

