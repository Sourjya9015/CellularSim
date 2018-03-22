classdef PLMod < hgsetget
    % PLMod:  Base class for a path loss model.
    % 
    % The path loss model is either between channel interfaces of
    % two different types (e.g. base stations and mobiles), or 
    % between channel interfaces of the same type (e.g. peer-to-peer
    % devices).  
    %
    % The base class stores a vector of channel interfaces for each type
    % using the map structure.  Channel interfaces should be added via the
    % addChanIf method.
    properties (Access = private)
        
        % Vector of channel interfaces and their names
        nTypes;
        chanIfNames;
        chanIfVecs;       
    end      
    
    properties (Constant = true, GetAccess = private)
        nTypesMax = 2;
    end
    
    methods
        
        % Constructor
        function obj = PLMod()  
            obj.nTypes = 0;
            obj.chanIfNames = cell(obj.nTypesMax,1);
            obj.chanIfVecs = cell(obj.nTypesMax,1);
        end
        
        % Destructor.  Clear references
        function delete(obj)
            for k = 1:obj.nTypesMax
                obj.chanIfVecs{k} = [];
            end
        end
        
        % Add channel interfaces
        function addChanIf(obj, name, chanIfVec)
            obj.nTypes = obj.nTypes + 1;
            if (obj.nTypes > obj.nTypesMax)
                error('Can only have %d channel interfaces for a path loss model', ...
                    obj.nTypesMax);
            end
            obj.chanIfNames{obj.nTypes} = name;
            obj.chanIfVecs{obj.nTypes} = chanIfVec;
        end
             
        % Get channel interface by name
        function chanIf = getChanIfByName(obj, name)
            exist = false;
            for ind = 1:obj.nTypes
                if (strcmp(name, obj.chanIfNames{ind}))
                    chanIf = obj.chanIfVecs{ind};
                    exist = true;
                end
            end
            if ~exist
                error('%s is not an interface type for the path loss model', ...
                    name);
            end
            
        end
        
        % Get channel interface by index
        function chanIf = getChanIfByInd(obj, ind)
            chanIf = obj.chanIfVecs{ind};
        end          
        
        % Accessors
        function nTypes = getNTypes(obj)
            nTypes = obj.nTypes;
        end
        function names = getNames(obj)
            names = obj.chanIfNames;
        end
        
        % Get index of a single channel element type        
        function [exist, ind] = getIndex(obj, name)
              exist = false;
            for ind = 1:obj.nTypes
                if (strcmp(name, obj.chanIfNames{ind}))
                    exist = true;
                    break;
                end
            end
        end
        
        % Get the order between two channel i/f types of names
        % name1 and name2.  If the class has only one i/f type, name1
        % must equal name2 and both must equal the name of the channel i/f
        % type.  If the class has two i/f types, name1 and name2 must be
        % the two i/f type name.  In this case, the output rev is true if
        % name1 and name2 are in "reverse" order to the way they are stored
        % in the map.
        function [exist,rev] = getOrder(obj, name1, name2)            
            rev = false;
            names = obj.chanIfNames;
            if (obj.nTypes == 1)
                exist = (strcmp(name1, names{1}) && strcmp(name2,names{2}));
            elseif (obj.nTypes ==2)
                if (strcmp(name1,names{1}) && strcmp(name2,names{2}))
                    exist = true;
                elseif (strcmp(name1,names{2}) && strcmp(name2,names{1}))
                    rev = true;
                    exist = true;
                else
                    exist = false;
                end
            else
                error('Path loss model must have one or two interface types.');
            end
        end
      
    end
    
    methods (Abstract)
        
        % Initialization function.  Called at the beginning to determine
        % shadowing etc.
        init(obj)        
        
        % Get path loss in dB between any sets of interfaces.
        % If ind1 and ind2 are unspecified or empty, then the method
        % returns a matrix pathLoss, where pl(i,j) is the path loss
        % between the i-th element of type name2 and j-th element of type
        % name 1.  The arguments ind1 and ind2 are subindices.
        pathLoss = getPathLoss(obj, name1, name2, ind1, ind2)
    end
    
end

