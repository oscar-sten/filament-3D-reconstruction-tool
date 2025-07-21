function G_minimal = removeBridgeNodes(G_minimal)

allDeg = degree(G_minimal);
nodeNames = G_minimal.Nodes.Name;

for i=1:numel(allDeg)

    % Identify the node's neighbors
    nghbrs = neighbors(G_minimal, nodeNames(i));
    % The node is a bridge node, if and only if it has two distinct
    % neighbors and has the degree 2. These criteria are not equivalent,
    % since a node can have multiple edges connecting to the same node
    % (i.e., the graph is not neccesarily a simple graph).
    if numel(nghbrs)==2 && degree(G_minimal, nodeNames(i)) == 2
        % Extract the weights of the bridge node's two edges.
        dists = distances(G_minimal, nodeNames(i), nghbrs);
        % Add a new edge between the two neighbors, the weight is the sum
        % of the two previous.
        G_minimal = addedge(G_minimal, nghbrs(1), nghbrs(2), (dists(1)+dists(2)));
        % Remove the bridge node from G_new
        G_minimal = rmnode(G_minimal, nodeNames(i));
    end

end



end