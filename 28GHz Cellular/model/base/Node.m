classdef Node < hgsetget
    % Node:  General class for a node
    properties (Access = private)
                       
        % Map of network devices
        % The key is the device name (a string, eg. 'lte-dev'),
        % The value is the NetDevice class.
        netDevMap = containers.Map();
        
        % Mobility model of class MobModel for determining the location of
        % the node.  Empty indicates that the position model is not set.
        mobMod = [];
    end
    
    methods
        
        % Constructor
        function obj = Node()            
        end
        
        % Destructor.  Nulls references for reference counting.
        function delete(obj)
            obj.mobMod = [];
            keys = obj.netDevMap.keys();
            for k = 1:length(keys)
                obj.netDevMap(keys{k}) = [];                
            end
        end
                   
        % Add a network device
        function installNetDev(obj, devName, dev)
            obj.netDevMap(devName) = dev;
            dev.install( obj );   
        end
        
        % Get the network device by name
        function dev = getNetDevByName(obj, devName)
            dev = obj.netDevMap(devName);
        end
        
        % Get network device by index.
        % If ind is not specified, it will get the first device.
        function dev = getNetDevByInd(obj, ind)
            if (nargin < 2)
                ind = 1;
            end
            k    = obj.netDevMap.keys();
            dev  = obj.netDevMap(k{ind});
        end
        
        % Set the mobility model
        function installMobMod(obj, mobMod)
            obj.mobMod = mobMod;
            obj.mobMod.install( obj );
        end
        
        % Get position model
        function mobMod = getMobMod(obj)
            mobMod = obj.mobMod;
        end             
     
    end
    
end