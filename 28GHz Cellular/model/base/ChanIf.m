classdef ChanIf < hgsetget
    % ChanIf:  Base class for a channel interface
    %
    % A channel interface operates between a network device (NetDev) and a
    % channel (Chan).  Each network device and each channel have one or 
    % more channel interfaces.    
    properties (Access = private) 
        % Network device (class NetDevice) for the channel interface.
        % Empty indicates it is unbound
        dev = [];
        
        % Channel that interface connects to.  
        chan = [];
    end
    
    methods
        
        % Constructor
        function obj = ChanIf()
        end         
        
        % Set and get device
        function setDev(obj, dev)
            obj.dev = dev;
        end
        function dev = getDev(obj)
            dev = obj.dev;
        end
        
        % Set and get channel
        function setChan(obj, chan)
            obj.chan = chan;
        end
        function chan = getChan(obj)
            chan = obj.chan;
        end
        
        % Clears references
        function delete(obj)
            obj.dev = [];
            obj.chan = [];
        end

        
    end
  
    
end