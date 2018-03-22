classdef MobMod < hgsetget     
     
    % MobMod:  Node mobility model base class
    properties (Access = private)
        % Node on which model is associated.
        node = [];
    end 

    methods 
        % Constructor
        function obj = MobMod()                       
        end
        
        % Called when the mobility model is installed on a node
        function install(obj, node)
            obj.node = node;
        end
        
        % Get reference to node
        function node = getNode(obj)
            node = obj.node;
        end
        
         % Destructor.  Clears references to avoid circular reference
         % counting.
        function delete(obj)
            obj.node = [];
        end
        
            
    end
    
end

