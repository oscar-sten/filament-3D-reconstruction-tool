function pairs = match_filament_crossing(G, nodeID1, nodeID2)

% Computes the arm vectors of the two input neighboring nodes and
% formulates an optimal matching problem for matching the arm nodes.
%
% Input: a graph object G with the node properties "Name" and "Coordinates".
%       nodeID1 and nodeID2, two 3-deg nodes.
% Output: pairs: [2 x 2] containing the matched neighbouring nodes.

node1 = findnode(G, nodeID1);
node2 = findnode(G, nodeID2);
coords1 = G.Nodes.Coordinates(node1, :);
coords2 = G.Nodes.Coordinates(node2, :);

% Path between pair
P = shortestpath(G, nodeID1, nodeID2);

% Compute crossing arm vectors
nghbrsNode1 = neighbors(G, nodeID1);
nghbrsNode1 = setdiff(nghbrsNode1, P);
armVectors1 = zeros(2,2);
nghbrsNode2 = neighbors(G, nodeID2);
nghbrsNode2 = setdiff(nghbrsNode2, P);
armVectors2 = zeros(2,2);
for i = 1:2
    % Node 1
    nghbrNode = findnode(G, nghbrsNode1{i});
    nghbrCoords = G.Nodes.Coordinates(nghbrNode, :);
    armVector =  nghbrCoords - coords1;
    armVectors1(:, i) = armVector';
    % Node 2
    nghbrNode = findnode(G, nghbrsNode2{i});
    nghbrCoords = G.Nodes.Coordinates(nghbrNode, :);
    armVector =  nghbrCoords - coords2;
    armVectors2(:, i) = armVector';

end

% Compute 2 x 2 cost matrix
costMatrix = zeros(2, 2);
for i=1:2
    for j=1:2
        v1 = armVectors1(:, i);
        v2 = armVectors2(:, j);
        costMatrix(i, j) = 180-rad2deg(compute_angle(v1, v2));

    end
end

% Solve optimal matching problem
matchings = matchpairs(costMatrix, 180);

% Assign the pairs. First row first pair (i.e., the index from the node1
% neighbors with assignment 1 and the index from the node2 neighbors
% with assignment 1. Corresponding assignment for row 2.
pairs = [nghbrsNode1(matchings(:, 1)==1), nghbrsNode2(matchings(:, 2)==1);
    nghbrsNode1(matchings(:, 1)==2), nghbrsNode2(matchings(:, 2)==2)];
end


