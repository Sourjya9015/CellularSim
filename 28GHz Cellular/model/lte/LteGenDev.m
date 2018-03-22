classdef LteGenDev < NetDevice 
    % LteGenDev:  A generic OFDMA device with multiple spectrum channel
    % interfaces
    properties (Access = private)
        
        
    end
    
    methods
        % Constructor
        % 
        % numIf = number of spectrum channel interfaces
        function obj = LteGenDev(ifName)
            % Network device
            obj = obj@NetDevice();
            
            % Create interfaces
            if (nargin >= 1)
                obj.createChanIf(ifName);
            end                    
            
        end                
        
        % Create spectrum channel interfaces
        function createChanIf(obj, ifName)
            
            % Spectrum channel interface
            if ~iscell(ifName)
                ifName = {ifName};
            end
            nif = length(ifName);
            specChanIf(nif,1) = SpecChanIf();            
            for iif = 1:nif
                obj.addChanIf(ifName{iif}, specChanIf(iif));
            end
            
        end                
                
    end
    
end

