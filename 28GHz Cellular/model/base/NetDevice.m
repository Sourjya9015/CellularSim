classdef NetDevice < hgsetget 
    % NetDevice:  Base class for a network device.
    % 
    % Each NetDevice gets installed on a node and connects to channels via
    % one one or more channel interfaces (class ChanIf).
    
    properties %(Access = private)
        
        % Map of channel interfaces
        % Key is a string (e.g. 'mmW-if')
        % Value if the channel interface
        % Multiple interfaces are used to support multiflow
        chanIfMap;
        
        % Node interface  -- this will be the socket interface to the device 
        nodeIf;
        
        % Node on which device is installed.
        % Empty indicates that the device is not installed
        node = [];
    end
    
    methods 
        % Constructor
        function obj = NetDevice()
            obj.chanIfMap = containers.Map();
        end
        
        % Destructor.  Explicitly clears references due to MATLAB problems
        % in garbage collection.
        function delete(obj)
            obj.node = [];
            keys = obj.chanIfMap.keys();
            for k = 1:length(keys)
                obj.chanIfMap(keys{k}) = [];
                obj.chanIfMap.remove(keys{k});
            end
        end

        
        % Add a channel interface
        function addChanIf(obj, name, chanIf)
            obj.chanIfMap(name) = chanIf;
            chanIf.setDev( obj );   % Set device on channel interface
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
        
        % Get a single channel interface.  Only works when there is only
        % one interface for the net device
        function chanIf = getChanIf(obj)
            if (obj.chanIfMap.Count() ~= 1)
                error('getChanIf assumes netDev has only 1 channel interface');
            end
            k       = obj.chanIfMap.keys();            
            chanIf  = obj.chanIfMap(k{1});
        end
        
        % Get all channel interfaces
        function chanIfMap = getChanIfMap(obj)
            chanIfMap = obj.chanIfMap;
        end
        
        % Called when a device is installed 
        function install(obj, node)
            obj.node = node;
        end
        
        % Get reference to node
        function node = getNode(obj)
            node = obj.node;
        end
        
    end
     
    
end