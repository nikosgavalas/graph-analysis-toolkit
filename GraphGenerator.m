classdef GraphGenerator
    methods
        function ret = complete(obj, N)
        % Constructor
        A = ones(N, N) - eye(N);
        ret = Graph(graph(A));
        end

        function ret = path(obj, N)
        % Generate path graph of N nodes
        range = 1:N-1;
        ret = Graph(graph(range, range + 1));
        end

        function ret = ring(obj, N)
        % Generate ring graph of N nodes
        first = 1:N;
        second = first + 1;
        second(end) = 1;
        ret = Graph(graph(first, second));
        end

        function ret = star(obj, N)
        % Generate star graph of N nodes
        ret = Graph(graph(ones(N - 1, 1), (1:N-1) + 1));
        end

        function ret = wheel(obj, N)
        % Generate wheel graph of N nodes
        % !use only for N>3!
        r = (1:N-1)+2;
        r(end) = 2;
        ret = Graph(graph( ...
            [ones(1,N-1) (1:N-1)+1], ...
            [(1:N-1)+1   r]));
        end

        function ret = random(obj, N, L, seed)
        % Generate a random graph with no self loops
        % !NOT necessarily connected!
        rng(seed);
        G = graph(true(N), 'omitselfloops');
        p = randperm(numedges(G), L);
        ret = Graph(graph(G.Edges(p, :)));
        end

        function ret = ER(obj, N, p)
        % Generate an Erdos-Renyi graph
        ret = ERGraph(N, p);
        end
    end
end
