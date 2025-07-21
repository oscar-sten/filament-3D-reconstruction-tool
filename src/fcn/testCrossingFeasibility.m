function flag = testCrossingFeasibility(G, nodeID1, nodeID2)

% Simple sanity test for the feasibility of an overlap. Helps in case of
% "oval" sturctures (i.e., two overlaps close to eachother).
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

if numel(nghbrsNode1)==2 && numel(nghbrsNode2)==2
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

    V1 = coords1 - coords2;

    V2 = coords2 - coords1;

    angle11 = rad2deg(compute_angle(V1, armVectors1(:, 1)));
    angle12 = rad2deg(compute_angle(V1, armVectors1(:, 2)));
    angle21 = rad2deg(compute_angle(V2, armVectors2(:, 1)));
    angle22 = rad2deg(compute_angle(V2, armVectors2(:, 2)));


    bool11 = dot(V1, armVectors1(:, 1))<0;
    bool12 = dot(V1, armVectors1(:, 2))<0;
    bool21 = dot(V2, armVectors2(:, 1))<0;
    bool22 = dot(V2, armVectors2(:, 2))<0;

    if bool11 || bool12 || bool21 || bool22
        flag = 1;
    else
        flag = 0;
    end
else
    flag=1;
end

end