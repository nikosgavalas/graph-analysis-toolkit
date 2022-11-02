clear all; clc; %clf;

%
% Graph class
%
% define specific graph through its links
g = graph( ...
    [1 1 2 2 2 3 3 3 4], ...
    [2 5 3 5 6 4 5 6 6] ...
    );
% then use the Graph class
G = Graph(g);
G.A;
%G.plot;
[upperSubgraph, lowerSubgraph] = G.subgraphsA(3);
G.N;
G.L;
G.d;
G.B;
G.Q;
G.avgDegree;
G.A_(2);
G.N_(2);
G.W_(2);
G.cc;
G.numTriangles;
G.distances;
G.hopcount;
G.diameter;
G.avgHopcount;
G.avgBetweeness;
G.degreeAssortativity;
G.X;
G.Lambda;
G.Z;
G.M;
G.numSpanningTrees;
G.algebraicConnectivity;
G.isConnected;
G.spectralRadius;

%
% ERGraph class
%
er = ERGraph(10, 0.1);
er.A;

%
% Electrical Network class
%
A = [[0 1 0 1];
     [1 0 1 1];
     [0 1 0 1];
     [1 1 1 0]];
en = ElectricalNetwork(graph(A), ones(5, 1));
x = [2; 0; 0; -2];
en.v(x);
en.y(x);
en.Omega;
en.R
en.closeness;
en.bestSpreader;

%
% Sample graphs to copy-paste and quickly test hypotheses
%
GG = GraphGenerator();
G1 = Graph(graph( ...      % a simple graph
    [1 1 1 2 2 2 3 3 3 4], ...
    [2 5 6 3 5 6 4 5 6 6] ...
    ));
G2 = GG.complete(5);       % complete graph
G3 = GG.path(5);           % path graph (is also tree. trees are a subclass of bipartites)
G4 = Graph(graph( ...      % bipartite graph
    [1 1 2 2 3 4], ...
    [4 6 3 5 4 5] ...
    ));
G5 = GG.ring(5);           % ring graph
G6 = Graph(graph( ...      % regular graph d=3, also bipartite
    [1 1 1 2 2 3 3 4 5 5 6 7], ...
    [2 4 5 3 6 4 7 8 6 8 7 8] ...
    ));
G7 = GG.star(5);           % star graph, also tree
G8 = GG.random(10, 17, 1); % random graph - this one's (seed=1) is connected
G9 = Graph(graph( ...      % small windmill graph (W(3,2))
    [1 1 2 3 3 3 3 4 6], ...
    [2 3 3 4 5 6 7 5 7]   ...
    ));
Gs = [G1 G2 G3 G4 G5 G6 G7 G8 G9];
names = ["Simple" "Complete" "Path" "Bipartite" "Ring" "Regular" "Star" "Random" "Windmill"];
res = zeros(1, size(Gs, 2));
for i = 1:size(Gs, 2)
    disp(names(i));
    g = Gs(i);
    % do stuff with g and quickly perform tests
    % example:
    %M = g.M;
    %res(i) = M(2, 2) > 0;
    % =================
    
    % =================
end
all(res == 1)

%
% Electrical Network class - always check the B matrix
%
en = ElectricalNetwork(graph( ...
    [[0 1 0 1];
     [1 0 1 1];
     [0 1 0 1];
     [1 1 1 0]]), ones(5, 1));
x = [2; 0; 0; -2];

%
% R-model simulation
%
GG = GraphGenerator();
N = 10;         % change me
LMAX = N * (N - 1) / 2;
k = LMAX - 20;  % change me
res1 = zeros(2 * k, 1);
res2 = zeros(2 * k, 1);
edges = zeros(k, 2);
g = GG.complete(N);
% remove edges
for i = 1:k
    g = g.getGraph;
    % save first edge
    e = g.Edges(1, 1).EndNodes;
    src = e(1);
    dst = e(2);
    edges(i, 1) = src;
    edges(i, 2) = dst;
    % remove it
    g = Graph(rmedge(g, 1));
    % get metrics
    %==========================
    res1(i) = g.spectralRadius;
    res2(i) = (g.L * 2) / g.N;
    %==========================
end
% add them back
for i = 1:k
    g = g.getGraph;
    src = edges(i, 1);
    dst = edges(i, 2);
    g = Graph(addedge(g, src, dst, 1));
    % get metrics
    %==========================
    res1(k + i) = g.spectralRadius;
    res2(k + i) = (g.L * 2) / g.N;
    %==========================
end
%plot(res1);
%hold on;
%plot(res2);
%legend('res1','res2'); % remember to clf
