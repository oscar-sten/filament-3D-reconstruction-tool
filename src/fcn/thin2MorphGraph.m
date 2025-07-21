function G_updated2 = thin2MorphGraph(G_updated)

G_updated2 = G_updated;

% Identify nodes with degree 2 and replace them with edges.

deg2Nodes = find(degree(G_updated2)==2);
bridgeNodes = G_updated2.Nodes.Name(deg2Nodes);

for i=1:numel(bridgeNodes)
    currentNode = bridgeNodes(i);
    nghbrs = neighbors(G_updated2, currentNode);
    if numel(nghbrs)==2 && degree(G_updated2, currentNode) == 2
        [~, dist] = shortestpath(G_updated2, nghbrs{1},  nghbrs{2});
        G_updated2 = addedge(G_updated2, nghbrs{1},  nghbrs{2}, dist);
        G_updated2 = rmnode(G_updated2, currentNode);
        % A node can have degree 2 without having two neighbors if 
        % the two edges are to the same node. This only happens in 
        % case of detection artifacts. 
    end
end

end