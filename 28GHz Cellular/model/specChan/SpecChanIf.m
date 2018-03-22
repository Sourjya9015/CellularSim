classdef SpecChanIf < ChanIf
    % SpecChanIf:  Interface to a spectrum channel
    properties 
        
        % Flags indicating if device performs TX and/or RX functions
        % Note that for FDD systems, there will be 
        txEn=false; 
        rxEn=false;
        
        % Tx parameters
        txPowdBm = 23;   % TX power in dBm
        
        % RX parameters
        noiseFig = 5;    % Noise figure in dB  
        
        % Antenna pattern that is applied in addition to the gain from the
        % array.
        sectInd = 1;        % sector index
        antPat = [];
        
        % Number of antennas in the horiz / vertical dimensions
        % Assumed to be a 2D array with lambda/2 spacing.
        nantH = 1;          
        nantV = 1;
        
     end
    
    methods
        
        % Constructor
        function obj = SpecChanIf(rxDir, txDir)
            obj = obj@ChanIf();
            if (nargin >= 1)
                obj.txDir = txDir;
            end
            if (nargin >= 2)
                obj.rxDir = rxDir;
            end
        end
        
        % Access direction
        function isTx = isTx(obj)
            isTx = obj.txDir;
        end
        function isRx = isRx(obj)
            isRx = obj.rxDir;
        end
        function setTxDir(obj, txDir)
            obj.txDir = txDir;
        end
        function setRxDir(obj, rxDir)
            obj.rxDir = rxDir;
        end
        
        % Add and get antenna pattern
        function addAntPattern(obj,antPat)
            obj.antPat = antPat;
            obj.antPat.sectInd = obj.sectInd;
        end
        
        % Get antenna gains for a vector of angles
        function gaindB = antGain(obj,angle)            
            if isempty(obj.antPat)
                nangle = length(angle);
                gaindB = zeros(nangle,1);
            else
                gaindB = obj.antPat.gain(angle);
            end
        end
                    
     
    end
    
end

