classdef SpecChan < Chan
    % SpecChan:  An OFDM channel
    % One channel is used for each band.  For FDD systems, there will be
    % two SpecChan classes.  For TDD systems, there will be one.
    properties % Public properties
        bwMHz = 10;     % Bandwidth in MHz
        fcMHz = 2100;   % Carrier frequency in MHz
    end
    properties (Access = protected)
        
        % Cell array of path loss models
        plMods = [];
    end
    
    % Constants
    properties (Constant=true)
        kT = -174;  % Thermal noise in dBm / Hz
    end
    
    methods
        % Constructor
        function obj = SpecChan()
            % Base Chan class
            obj = obj@Chan();
        end
        
        % Add a path loss model
        function addPlMod(obj, plMod)
            obj.plMods = cat(1, obj.plMods, plMod);
        end   
        
        % Get path loss in dB between any sets of interfaces types.
        % If ind1 and ind2 are unspecified or empty, then the method
        % returns a matrix pathLoss, where pl(i,j) is the path loss
        % between the i-th element of type name2 and j-th element of type
        % name 1.  The arguments ind1 and ind2 are subindices.
        function pathLoss = getPathLoss(obj, name1, name2, ind1, ind2)       
            % Loop over path loss models
            found = false;
            nmod = length(obj.plMods);
            for imod = 1:nmod
                mod = obj.plMods(imod);
                found = mod.getOrder(name1,name2);
                if (found)
                    if (nargin < 4)
                        ind1 = [];
                    end
                    if (nargin < 5)
                        ind2 = [];
                    end
                    pathLoss = mod.getPathLoss(name1,name2,ind1,ind2);
                    break;
                end
            end
            if ~found
                error('No path loss models between %s and %s', name1, name2);
            end
        end
        
        function plMod = getPlMod(obj, name1, name2)       
            % Loop over path loss models
            found = false;
            nmod = length(obj.plMods);
            for imod = 1:nmod
                mod = obj.plMods(imod);
                found = mod.getOrder(name1,name2);
                if (found)
                    plMod = mod;
                    break;
                end
            end
            if ~found
                error('No path loss models between %s and %s', name1, name2);
            end
        end
        
        % Compute SINR assuming cell selection to the strongest cell
        function [sinr, Icell, pathLoss, txpow] = computeSinr(obj, rxname, txnames)
            % Loop over all the TX types
            if ~iscell(txnames)
                txnames = {txnames};
            end
            ntxTypes = length(txnames);
            for itx = 1:ntxTypes                     
                txnamei = txnames{itx};                
                pathLossi = obj.getPathLoss(txnamei,rxname);
                chanIf = obj.getChanIfByName(txnamei);
                n = length(chanIf);
                powi = zeros(n,1);
                for i = 1:n
                    powi(i) = chanIf{i}.txPowdBm;
                end
                if (itx == 1)
                    txpow = powi;
                    pathLoss = pathLossi;
                else
                    txpow = [txpow; powi];
                    pathLoss = [pathLoss pathLossi];
                end
            end
            
            % Get noise power
            chanIf = obj.getChanIfByName(rxname);
            nrx = size(pathLoss,1);
            noisePow = zeros(nrx,1);
            ktPow = obj.kT + 10*log10(obj.bwMHz) + 60;
            for irx = 1:nrx
                noisePow(irx) = 10^(0.1*(chanIf{irx}.noiseFig + ktPow));
            end
            
            % Compute SINR
            p = 10.^(0.1*(repmat(txpow',nrx,1) - pathLoss));
            [sigPow,Icell] = max(p,[],2);   % nb: Icell/Isel: max p !
            totPow = sum(p,2) + noisePow;
            sinr = sigPow./(totPow - sigPow);            
            sinr = 10*log10(sinr);
            
        end
        
    end
    
end

