classdef Chan < hgsetget
    % Chan:  Base class for a channel
    %
    % A channel represents the medium by through which multiple network
    % devices communicate.  The Chan class may instatiate multiple
    % point-to-point channels -- it is not intended to be a single link.
    properties (Access = private)
        
        % Channel interface vector map.
        % The "key" is the type of channel interface,
        % eg. 'TX' and 'RX' or 'eNB' and 'UE'.
        % Each "value" is a vector of channel interfaces for this interface
        % type.
        chanIfMap = containers.Map();
        
    end
    
    methods
        
        % Constructor
        function obj = Chan()            
        end
        
        % Destructor.  Clear references
        function delete(obj)
            names = obj.chanIfMap.keys();
            for k = 1:length(names)
                obj.chanIfMap(names{k}) = [];
            end
        end
        
        % Add channel interfaces
        function addChanIf(obj, name, chanIfVec)
            obj.chanIfMap(name) = chanIfVec;
        end
             
        % Get channel interface by name
        function chanIf = getChanIfByName(obj, name)
            chanIf = obj.chanIfMap(name);
        end
        
        % Get channel interface by index
        function chanIf = getChanIfByInd(obj, ind)
            if (nargin < 2)
                ind = 1;
            end
            k       = obj.chanIfMap.keys();
            chanIf  = obj.chanIfMap(k{ind});
        end
        
        % Get all channel interface vectors
        function chanIfMap = getChanIfMap(obj)
            chanIfMap = obj.chanIfMap;
        end

             
    end
    
end