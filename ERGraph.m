classdef ERGraph < Graph
    %ERGRAPH Erdos-Renyi graph
    %   Inherits the base Graph class
    properties
        p
    end

    methods
        function obj = ERGraph(N, p)
        % Constructor
        A = rand(N, N) > p;
        A = not(triu(ones(N, N))) .* A;
        A = A + A';
        obj@Graph(graph(A));
        obj.p = p;
        end
    end
end
