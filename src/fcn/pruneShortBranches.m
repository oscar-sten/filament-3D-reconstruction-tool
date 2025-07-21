function G = pruneShortBranches(SG)

G_minimal = graph(string(SG.Edges.EndNodes(:, 1)),...
    string(SG.Edges.EndNodes(:, 2)),...
    SG.Edges.Weight);

% Get the degree of each node
nodeDegrees = degree(G_minimal);

% Find nodes with degree 1
nodesDegree1 = find(nodeDegrees == 1);

nodesDegree1_ID = G_minimal.Nodes.Name(nodesDegree1);

% For each node with degree 1, check if its neighbor has degree 3
for i = 1:length(nodesDegree1)
    node = nodesDegree1_ID(i);
    nghbrs = neighbors(G_minimal, node);
    
    % If the neighbor has degree 3, remove the node
    if numel(nghbrs) == 1 && degree(G_minimal, nghbrs) == 3
        G_minimal = rmnode(G_minimal, node);
    end
end

G = G_minimal;

end