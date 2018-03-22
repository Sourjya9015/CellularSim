classdef HetNet < hgsetget
    % HetNet:  Base class for a heterogenous network
    properties
        
        % Data on the nodes and their channel interfaces
        nodeDatMap;
        chanDatMap;
        chanIfDatMap;
        
        % Maps of the various elements in the network
        nodeMap;        % nodes, indexed by node names
        mobModMap;      % mobility models, indexed by node names
        devMap;         % devices, indexed by chan i/f names
        chanIfMap;      % channel interfaces
        chanMap;        % channels, indexed by chan names
        plModMap;       % path loss models
               
        % Device names for different sectors
        devNameSect;
        
        % Position limit
        posLim = [2000 2000]';                 
        
    end
    
    methods
        
        % Constructor.  Note that initialization is done in init so that
        % the properties can be set first
        function obj = HetNet(nodeDatMap, chanDatMap, chanIfDatMap)
            
            % Get names and number of each type
            nodeNames = nodeDatMap.keys();
            chanNames = chanDatMap.keys();
            chanIfNames = chanIfDatMap.keys();
            ntNode = length(nodeNames);
            ntChan = length(chanNames);
            ntChanIf = length(chanIfNames);
            
            % Save the data maps
            obj.nodeDatMap = nodeDatMap;
            obj.chanDatMap = chanDatMap;
            obj.chanIfDatMap = chanIfDatMap;
                        
            % Create empty maps     
            obj.plModMap = containers.Map();
            obj.nodeMap = containers.Map();
            obj.devMap = containers.Map();
            obj.chanIfMap = containers.Map();
            obj.mobModMap = containers.Map();
            obj.chanMap = containers.Map();
                        
            % Create the channels                       
            for it = 1:ntChan
                obj.chanMap(chanNames{it}) = SpecChan();
            end
            
            % Create the nodes
            for it = 1:ntNode
                
                % Create the nodes and mobility models
                nodeName = nodeNames{it};
                nodeDat = obj.nodeDatMap(nodeName);
                num = nodeDat.num;
                nodeList = cell(num,1);
                mobModList = cell(num,1);
               
                % Install the mobility models in the nodes
                for inode = 1:num  
                    nodeList{inode} = Node();
                    mobModList{inode} = ConstMobMod();
                    nodeList{inode}.installMobMod( mobModList{inode} );
                end
                obj.nodeMap(nodeName) = nodeList;
                obj.mobModMap(nodeName) = mobModList; 
     
            end
            
            % Create the devices and channel interfaces
            for it = 1:ntChanIf
                
                chanIfName = chanIfNames{it};
                chanIfDat = obj.chanIfDatMap(chanIfName);
                nodeName = chanIfDat.nodeName;
                chanName = chanIfDat.chanName;
                nnode = nodeDatMap(nodeName).num;
                nodeList = obj.nodeMap(nodeName);
                
                % Check if it is sectorized
                if (isfield(chanIfDat,'nsect'))
                    nsect = chanIfDat.nsect;
                    ifName = cell(nsect,1);                    
                    for isect = 1:nsect
                        ifName{isect} = sprintf('%s:%d', chanName, isect);
                    end
                else
                    nsect = 1;
                    ifName = {chanName};
                end
                
                % Create the devices and channel interfaces
                ndev = nsect*nnode;
                devList = cell(ndev,1);
                chanIfList = cell(ndev,1);
       
                % Install the devices in the nodes
                idev = 0;
                for inode = 1:nnode
                    for isect = 1:nsect
                        idev = idev+1;
                        dev = LteGenDev();
                        nodeList{inode}.installNetDev( ifName{isect}, dev);
                        dev.createChanIf( chanIfName );
                        devList{idev} = dev;
                        chanIf = dev.getChanIf();
                        chanIf.sectInd = isect;
                        chanIfList{idev} = chanIf;
                        
                    end
                end

                % Save devices and channel interfaces
                obj.devMap(chanIfName) = devList;
                obj.chanIfMap(chanIfName) = chanIfList;
                
                % Add the interfaces to the channel
                chan = obj.chanMap(chanName);
                chan.addChanIf( chanIfName, chanIfList );
                               
            end
                 
        end
        
        % Destructor.  Clear references in the maps
        function delete(obj)                       
            maps = {obj.plModMap, obj.nodeMap, obj.devMap, ...
                obj.chanIfMap, obj.chanMap, obj.plModMap};
            for imap = 1:length(maps)
                m = maps{imap};
                names = m.keys();
                for k = 1:length(names)
                    m(names{k}) = [];
                end
            end                    
            
        end              
        
        % Create a random placement.
        function drop(obj)
    
            % Randomly drop the nodes
            nodeNames = obj.nodeMap.keys();
            ntNode = length(nodeNames);
            for it = 1:ntNode
                
                % Get the mobility model
                nodeName = nodeNames{it};
                mobModList = obj.mobModMap(nodeName);
                
                % Get the reference to layout model
                layoutMod = obj.nodeDatMap(nodeName).layoutMod;
                nnode = length(mobModList);
                loc = layoutMod.drop(nnode);
                
                % Set the locations from the             
                for inode = 1:nnode                   
                    mobModList{inode}.setLoc( loc(inode,:) );
                end
            end
            
            % Initialize all the path loss models
            modNames = obj.plModMap.keys();
            nmod = length(modNames);
            for imod = 1:nmod
                plMod = obj.plModMap(modNames{imod});
                plMod.init();
            end
        end
        
        % Get channel by name.  If name is not specified, it gets the first
        % channel
        function chan = getChan(obj, chanName)
            if (nargin < 2)
                useFirst = true;
            else
                useFirst = isempty(chanName);
            end
            if (~useFirst)
                chan = obj.chanMap(chanName);            
            else
                k = obj.chanMap.keys();
                chan = obj.chanMap(k{1});
            end
        end
        
        % Get list of nodes by node name
        function nodeList = getNodeList(obj, nodeName) 
            nodeList = obj.nodeMap(nodeName);
        end
        % Get list of mobility models by node name
        function mobModList = getMobModList(obj, nodeName) 
            mobModList = obj.mobModMap(nodeName);
        end

        % Get list of devices by channel interface name
        function devList = getDevList(obj, chanIfName) 
            devList = obj.devMap(chanIfName);
        end
        % Get list of channel interfaces by channel interface name
        function chanIfList = getChanIfList(obj, chanIfName) 
            chanIfList = obj.chanIfMap(chanIfName);
        end
                
        % Add a path loss model
        function addPlMod(obj, plModName, chanIfName1, chanIfName2, opt)
            
            plModTypeflag = false;
            chanIf2flag = false;
            if nargin>4
                plModTypeflag = true;
                chanIf2flag = true;
            end
            
            if nargin==4
                if isfield(chanIfName2,'plModType')
                    plModTypeflag = true;
                else
                    chanIf2flag = true;
                end
            end
            
            plModType = 'dist';
            if plModTypeflag
                plModType = opt.plModType;
            end
            % Create a distance-based path loss model
            switch plModType
                case 'dist'
                    plMod = DistPlMod();
                case 'hybrid'
                    plMod = HybDistPlMod();
                case 'mimo'
                    plMod = MimoPlMod();
            end
            
            % Add the channel interfaces
            chanIf1 = obj.getChanIfList(chanIfName1);
            plMod.addChanIf(chanIfName1, chanIf1);
            if chanIf2flag
                chanIf2 = obj.getChanIfList(chanIfName2);
                plMod.addChanIf(chanIfName2, chanIf2);
            end
            
            % Set the position limits
            plMod.setSqWrap(obj.posLim(1), obj.posLim(2));
            
            % Add the path loss model
            obj.plModMap(plModName) = plMod;     
            
            % Add the path loss model to the channel
            chanIfDat = obj.chanIfDatMap(chanIfName1);
            chanName = chanIfDat.chanName;
            chan = obj.getChan(chanName);
            chan.addPlMod(plMod);
                
        end
        
        % Get path loss model by name
        function plMod = getPlMod(obj, plModName)
            plMod = obj.plModMap(plModName);
        end
            
    end
    
end
