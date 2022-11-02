classdef Graph
    %GRAPH Wrapper class for MATLAB graphs
    %   Contains shortcuts for graph metrics

    properties
        graph
    end

    methods
        function obj = Graph(graph)
        % Constructor
        obj.graph = graph;
        end

        function ret = N(obj)
        % Number of nodes
        ret = numnodes(obj.graph);
        end

        function ret = L(obj)
        % Number of links (edges)
        ret = numedges(obj.graph);
        end

        function ret = A(obj)
        % Adjacency matrix
        ret = full(adjacency(obj.graph));
        end

        function ret = plot(obj)
        % Plot the graph
        ret = plot(obj.graph);
        end

        function [upper, lower] = subgraphsA(obj, splitPoint)
        % Split into two subgraphs
        adj = obj.A;
        upper = adj(1:splitPoint-1, 1:splitPoint-1);
        lower = adj(splitPoint:obj.N, splitPoint:obj.N);
        end

        function ret = d(obj)
        % Degree vector
        ret = degree(obj.graph);
        end

        function ret = degree(obj)
        % Degree vector alias
        ret = obj.d;
        end

        function ret = B(obj)
        % Incidence matrix
        ret = full(incidence(obj.graph));
        end

        function ret = Q(obj)
        % Laplacian
        ret = full(laplacian(obj.graph));
        end

        function ret = avgDegree(obj)
        % Average node degree
        ret = 2 * obj.L / obj.N;
        end

        function ret = A_(obj, k)
        % Number of k-hop walks between all nodes
        ret = obj.A ^ k;
        end
        
        function ret = N_(obj, k)
        % Total number of k-hop walks
        ret = sum(sum(obj.A_(k)));
        end

        function ret = W_(obj, k)
        % Total number of closed k-hop walks in G
        ret = trace(obj.A_(k));
        end

        function ret = cc(obj)
        % Clustering coefficient
        d = obj.d;
        ret = obj.W_(3) / (d' * d - 2 * obj.L);
        end

        function ret = numTriangles(obj)
        % Number of triangles
        ret = (1 / 6) * (obj.N_(2) - obj.W_(2)) * obj.cc;
        end

        function ret = distances(obj)
        % Shortest distances between nodes
        ret = distances(obj.graph);
        end

        function ret = hopcount(obj)
        % Hopcount = path length between nodes.
        % Essentially alias of distances
        ret = obj.distances;
        end

        function ret = diameter(obj)
        % Diameter
        ret = max(max(obj.distances));
        end

        function ret = diameter_(obj, limit)
        % Do not use. Use diameter() instead
        % Returns -1 if hopcount > limit
        ret = -1;
        A = obj.A;
        I = eye(obj.N);
        for i = 1:limit
            if all(all((I + A)^i > 0))
                ret = i;
                return;
            end
        end
        end

        function ret = avgHopcount(obj)
        % Average hopcount
        ret = sum(sum(triu(obj.distances))) / ((obj.N - 1) * obj.N / 2);
        end

        function ret = avgBetweeness(obj)
        % Average betweeness
        ret = (1 / obj.L) * (obj.N * (obj.N - 1) / 2) * obj.avgHopcount;
        end

        function ret = degreeAssortativity(obj)
        % Degree assortativity
        N1 = obj.N_(1);
        N2 = obj.N_(2);
        N3 = obj.N_(3);
        sd3 = sum(obj.d .^ 3);
        ret = (N1 * N3 - N2^2) / (N1 * sd3 - N2^2);
        end

        function ret = Lambda(obj)
        % Eigenvalues of the adj matrix
        [~, Lambda] = eig(obj.A);
        ret = Lambda;
        end

        function ret = X(obj)
        % Eigenvectors of the adj matrix
        [X, ~] = eig(obj.A);
        ret = X;
        end

        function ret = Z(obj)
        % Eigenvectors of the laplacian
        [Z, ~] = eig(obj.Q);
        ret = Z;
        end

        function ret = M(obj)
        % Eigenvalues of the laplacian
        [~, M] = eig(obj.Q);
        ret = M;
        end

        function ret = numSpanningTrees(obj)
        % Number of spanning trees
        M = obj.M;
        ret = (1 / obj.N) * prod(diag(M(2:end, 2:end)));
        end

        function ret = algebraicConnectivity(obj)
        % Algebraic connectivity
        M = obj.M;
        ret = M(2, 2);
        end

        function ret = isConnected(obj)
        % returns one if the graph is connected
        ret = all(conncomp(obj.graph) == 1); % the matlab way
        %ret = obj.algebraicConnectivity > 10e-10; % > 0 normally but due to
        %numerical errors I should use something like 10e-10
        end

        function ret = closeness(obj)
        % Closeness matrix. ret_i is the closeness Cl_i of the i-th node.
        ret = 1 ./ sum(obj.distances);
        end

        function ret = spectralRadius(obj)
        % Spectral radius
        Lambda = obj.Lambda;
        ret = Lambda(obj.N, obj.N);
        end

        function ret = getGraph(obj)
        % Return the underlying graph
        ret = obj.graph;
        end
    end
end
