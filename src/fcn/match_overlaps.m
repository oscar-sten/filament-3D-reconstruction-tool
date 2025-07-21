function G_updated = match_overlaps(G_input, G_minimal,...
    nodeIDCrossingPairs)
% Returns the updated representative/thin graph with overlaps.
% Input: 
%   G_input: current representative/thin graph.
%   G_minimal: the minimal homeomorphic graph without overlaps.
%   nodeIDCrossingPairs: [N x 2], ID:s of branch nodes determined to be
%   part of an overlap structure.

G_updated = G_input;

n_crossing_pairs = size(nodeIDCrossingPairs, 1);

nghbrs2Connect = [];

for i=1:n_crossing_pairs
    nodeID1 = string(nodeIDCrossingPairs(i, 1));
    nodeID2 = string(nodeIDCrossingPairs(i, 2));
    pairs = match_filament_crossing(G_input, nodeID1, nodeID2);

    % Candidates
    pair1 = pairs(1, :);
    pair2 = pairs(2, :);
    
    middleRegion = shortestpath(G_input, nodeID1, nodeID2);
    middleRegion_idx = findnode(G_input, middleRegion);
    % Best option
    middleRegion_depth = G_input.Nodes.Z_coordinates(middleRegion_idx);
    
    % Second best option
    idx_invalid = middleRegion_depth==-1;
    if sum(idx_invalid)==0
        middleRegion_depth(idx_invalid) = G_input.Nodes.SourceGraphIdx(...
            middleRegion_idx(idx_invalid));
    else
        middleRegion_depth = middleRegion_depth(~idx_invalid);
    end

    % Determine depth of candidates
    nodePair1 = findnode(G_updated, pair1);
    nodePair2 = findnode(G_updated, pair2);
    
    % Best option
    depthPair1 = G_updated.Nodes.Z_coordinates(nodePair1);
    depthPair2 = G_updated.Nodes.Z_coordinates(nodePair2);

    % Second best option
    idx_invalid = depthPair1==-1;
    if sum(idx_invalid)==0
        depthPair1(idx_invalid) = G_input.Nodes.SourceGraphIdx(...
            nodePair1(idx_invalid));
    else
        depthPair1 = depthPair1(~idx_invalid);
        
    end

    idx_invalid = depthPair2==-1;
    if sum(idx_invalid)==0
    depthPair2(idx_invalid) = G_input.Nodes.SourceGraphIdx(...
        nodePair2(idx_invalid));
    else
        depthPair2 = depthPair2(~idx_invalid);
    end

    pair1_diff = abs(mean(depthPair1)-mean(middleRegion_depth));
    pair2_diff = abs(mean(depthPair2)-mean(middleRegion_depth));

    if mean(pair1_diff) < mean(pair2_diff)
        nghbrs2Connect = [nghbrs2Connect; pair2];
    else
        nghbrs2Connect = [nghbrs2Connect; pair1];
    end

end


traversals = traversalsFromPairs(nodeIDCrossingPairs);

% Merger if needed
for i=1:numel(traversals)
    if numel(traversals{i})>2
        current_traversal = traversals{i};
        node1 = findnode(G_input, string(current_traversal(1)));
        node2 = findnode(G_input, string(current_traversal(end)));
        [P, ~] = shortestpath(G_input, node1, node2);
        mergers = intersect(G_input.Nodes.Name(P), nghbrs2Connect(:));
        nodes2Merge = zeros(1, 2);
        rows2remove = [];
        for j = 1:numel(mergers)
            [x, y] = find(nghbrs2Connect==mergers(j));
            not_y = setdiff([1 2], y);
            corresp = nghbrs2Connect(x, not_y);
            nodes2Merge(1, j) = corresp;
            rows2remove = [rows2remove; x];
        end

        if ~any(nodes2Merge==0)
            nghbrs2Connect = [nghbrs2Connect; nodes2Merge];
            nghbrs2Connect(rows2remove, :) = [];
        end
    end
end


% Merge neighbors

n_pairs2Merge = size(nghbrs2Connect, 1);

for i=1:n_pairs2Merge
    % Identify nodes and neighbors
    nodeID1 = nghbrs2Connect(i, 1);
    nodeID2 = nghbrs2Connect(i, 2);
    nghbrs1 = neighbors(G_updated, nodeID1);
    nghbrs2 = neighbors(G_updated, nodeID2);
    % Identify endnodes of edges to remove
    [P, weight] = shortestpath(G_input, nodeID1, nodeID2);
    bp1 = intersect(nghbrs1, P);
    bp2 = intersect(nghbrs2, P);

    % Remove edge to BP
    G_updated = rmedge(G_updated, nodeID1, bp1);
    G_updated = rmedge(G_updated, nodeID2, bp2);

    % Add edges to merge
    G_updated = addedge(G_updated, nodeID1, nodeID2, weight);

end


% Process 4-connected neighbors

degrees = degree(G_minimal);
nodes_degree_4 = find(degrees == 4); % nodes of degree 4
nodes_degree_4IDs = G_minimal.Nodes.Name(nodes_degree_4);
perfect_crossing_coords = G_minimal.Nodes.Coordinates(nodes_degree_4, :);

n_4connectedNodes = numel(nodes_degree_4IDs);
for k = 1:n_4connectedNodes

    centerNode = findnode(G_updated, nodes_degree_4IDs(k));

    fourNodeNghbrs = neighbors(G_updated, centerNode);

    if numel(fourNodeNghbrs)==4 % Check that the 4 neighbors are actually 4

        centerCoords = G_updated.Nodes.Coordinates(centerNode, :);

        nghbrCoords = G_updated.Nodes.Coordinates(fourNodeNghbrs, :);

        armVectors = zeros(4, 2);

        for i=1:4
            armVectors(i, :) = centerCoords - nghbrCoords(i, :);
        end

        costMatrix = zeros(4, 4);

        for i = 1:4
            for j=1:4
                costMatrix(i, j) = 180-rad2deg(compute_angle(armVectors(i, :),...
                    armVectors(j, :)));
            end
        end

        pairs = matchpairs(costMatrix, 180);

        pair2merge = fourNodeNghbrs(pairs(1, :));

        pair2mergeIDs = G_updated.Nodes.Name(pair2merge);

        centerNodeID = G_updated.Nodes.Name(centerNode);

        G_updated = rmedge(G_updated, centerNodeID, pair2mergeIDs(1));
        G_updated = rmedge(G_updated, centerNodeID, pair2mergeIDs(2));

        % Add edges to merge
        [~, weight] = shortestpath(G_updated, pair2mergeIDs(1), pair2mergeIDs(2));
        G_updated = addedge(G_updated, pair2mergeIDs(1), pair2mergeIDs(2), weight);

    end
end

end

