classdef PLModFix < PLMod
    % PLFixMod:  Path loss model based on fixed path loss    
    properties (Access = private)
                  
        % Path loss matrix plMat(itx,irx) = path loss in dB from TX itx to
        % RX irx
        plMat;  
                
    end
    
    methods
        
        % Constructor
        function obj = PLModFix(ifList, plMat)
            
            % Add the network devices
            obj = obj@PLMod(ifList);                 
            
            % Save the path loss matrix
            [nrx,ntx] = size(plMat);
            [Irx,Itx] = obj.getInd();            
            if ((nrx ~= length(Irx)) || (ntx ~= length(Itx)))
                error('Path loss matrix is not correct size');
            end
            obj.plMat = plMat;
        end           
        
        % Get path loss matrix
        function plMat = getPlMat(obj)
            plMat = obj.plMat;
        end
        
        % Path loss function
        function pl = getPathLoss(obj, irx, itx)
            pl = obj.plMat(irx,itx);
        end
        
    end
    
end

