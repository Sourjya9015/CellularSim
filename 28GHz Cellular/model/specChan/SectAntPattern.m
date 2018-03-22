classdef SectAntPattern < AntPattern
    % SectAntPattern:  Antenna pattern as specfied in 36.814    
    properties
        secWid = 2*pi*7/36; % sector width in angles
        angle0 = 0;      % bore sight angle for sectInd = 1
  
        % Pattern type
        patType = '3GPP';
        
        % Parameters for '3GPP' model (36.814)
        % gain = -min(12*(angle/angle3dB)^2,Amin)        
        angle3dB = 70*pi/180;
        Amin = 25;
        
        % Parameters for 'flat' model:
        % gain = gainMax        in sector
        %      = gainMax-fbGain outside sector
        fbGain = 80;    % front-to-back gain
        gainMax = 0;    % max gain      
        
    end
    
    methods
        % Constructor
        function obj = SectAntPattern(nsect,fbGain,gainMax,angle0)
            obj = obj@AntPattern();
            if (nargin >= 1)
                obj.secWid = 2*pi/nsect;
            end
            if (nargin >= 2)
                obj.fbGain = fbGain;
            end
            if (nargin >= 3)
                obj.gainMax = gainMax;
            end
            if (nargin >= 4)
                obj.angle0 = angle0;
            end
        end
        
        % Returns gain as a function of the angle in radians.
        % Note that angle may be a vector of angles
        function gaindB = gain(obj, angle)
            dangle = mod(angle - obj.sectInd*obj.secWid - obj.angle0, 2*pi);
            if (strcmp(obj.patType, '3GPP'))
                gaindB = -min(12*(dangle/obj.angle3dB).^2, obj.Amin);
            elseif (strcmp(obj.patType, 'flat'))           
                I = ~((dangle < obj.secWid/2) | (dangle > 2*pi-obj.secWid/2));
                gaindB = obj.gainMax - I*obj.gainMax - I*obj.fbGain;
            else
                error('Unknown antenna pattern type %s', obj.patType);                
            end
        end
        
    end
    
end

