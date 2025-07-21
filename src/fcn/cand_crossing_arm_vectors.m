function [armVectors, neighborNodes] = cand_crossing_arm_vectors(...
    G_minimal, G_thin, nodePair)


node1 = nodePair(1);
node2 = nodePair(2);
coords1 = G_minimal.Nodes.Coordinates(node1, :);
coords2 = G_minimal.Nodes.Coordinates(node2, :);
nodeID1 = G_minimal.Nodes.Name(node1);
nodeID2 = G_minimal.Nodes.Name(node2);


% Path between pair
P = shortestpath(G_thin, nodeID1, nodeID2);


nghbrsNode1 = neighbors(G_thin, nodeID1);
nghbrsNode1 = setdiff(nghbrsNode1, P);


nghbrsNode2 = neighbors(G_thin, nodeID2);
nghbrsNode2 = setdiff(nghbrsNode2, P);

armVectors = zeros(2, 2, 2);
neighborNodes = zeros(2, 2);

for j = 1:2
    % Node 1
    nghbrNode = findnode(G_thin, nghbrsNode1{j});
    nghbrCoords = G_thin.Nodes.Coordinates(nghbrNode, :);
    armVector = coords1 - nghbrCoords;
    armVectors(j, :, 1) = armVector;
    neighborNodes(j, 1) = string(nghbrsNode1{j});

    % Node 2
    nghbrNode = findnode(G_thin, nghbrsNode2{j});
    nghbrCoords = G_thin.Nodes.Coordinates(nghbrNode, :);
    armVector = coords2 - nghbrCoords;
    armVectors(j, :, 2) = armVector;
    neighborNodes(j, 2) = string(nghbrsNode2{j});
end


end