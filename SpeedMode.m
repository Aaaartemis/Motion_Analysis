classdef SpeedMode
    properties
        Mode
        Time
    end
    
    methods
        % Constructor Method
        function obj = SpeedMode(mode, time)
            obj.Mode = mode;
            obj.Time = time;
        end
        
    end
end
