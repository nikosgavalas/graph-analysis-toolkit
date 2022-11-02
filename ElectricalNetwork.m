classdef ElectricalNetwork < Graph
    %ELECTRICALNETWORK Electrical Network graph
    %   Inherits the base Graph class

    properties
        resistances
    end

    methods
        function obj = ElectricalNetwork(graph, resistances)
        % Constructor
        obj@Graph(graph);
        obj.resistances = resistances;
        end

        function ret = v(obj, x)
        % Voltages, x is the injected currents vector
        % the result returned assumes v_avg = 0. If you want the result
        % with reference to a specific node (so v of that node == 0), you
        % must add this potential to the whole vector. For example, if
        % v = [0.7500 0.0000 -0.2500 -0.5000] and you want v_4 == 0, add
        % -0.5 to the vector so you get [1.2500 0.5000 0.2500]
        ret = pinv(obj.Q) * x;
        end

        function ret = y(obj, x)
        % Current vector
        % x is the injected currents vector
        B = obj.B;
        ret = diag(1 ./ obj.resistances) * B' * obj.v(x);
        end

        function ret = Q(obj)
        % Weighted laplacian, overrides the normal laplacian
        B = obj.B;
        ret = B * diag(1 ./ obj.resistances) * B';
        end

        function ret = zeta(obj)
        % Diagonal of the pseudoinverse of the Laplacian
        ret = diag(pinv(obj.Q));
        end

        function ret = Omega(obj)
        % Effective graph resistance matrix
        u = ones(obj.N, 1);
        zeta = obj.zeta;
        ret = u * zeta' + zeta * u' - 2 * pinv(obj.Q);
        end

        function ret = R(obj)
        % Effective graph resistance R_G
        ret = obj.N * trace(pinv(obj.Q));
        % other ways:
        %ret = (1 / 2) * sum(sum(obj.Omega));
        %M = obj.M
        %ret = obj.N * sum(1 ./ diag(M(2:end, 2:end)));
        end

        function ret = bestSpreader(obj)
        % Best spreader node
        % better look at the whole zeta, because there may be more than one
        % nodes with equal values
        [~, argmin] = min(obj.zeta);
        ret = argmin;
        end
    end
end
