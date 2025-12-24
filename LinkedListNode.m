classdef LinkedListNode < handle
    properties
        Data
        Next
    end
    
    methods
        % Constructor Method
        function obj = LinkedListNode(data)
            obj.Data = data;
            obj.Next = LinkedListNode.empty;
        end
        
    end
end
