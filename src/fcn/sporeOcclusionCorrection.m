function G_updated = sporeOcclusionCorrection(G_updated, ...
    nodeSets, armVectorSet, connectionVectorSet)

for i=1:numel(nodeSets)
    currentNodeSet = nodeSets{i};
    if numel(currentNodeSet) > 1
        armVectors = armVectorSet{i};
        connectionVectors = connectionVectorSet{i};
        combos = nchoosek(1:numel(currentNodeSet), 2);
        n_combos = size(combos, 1);
        point_angles = zeros(n_combos, 2);
        for j=1:n_combos
            node1 = currentNodeSet(combos(j, 1));
            node2 = currentNodeSet(combos(j, 2));
            % Arm vectors
            armVector1 = armVectors(combos(j, 1), :);
            armVector2 = armVectors(combos(j, 2), :);
            % Connection vectors
            ordinal1 = setdiff(1:numel(currentNodeSet), combos(j, 1));
            ordinal1 = (combos(j, 2)==ordinal1);
            ordinal2 = setdiff(1:numel(currentNodeSet), combos(j, 2));
            ordinal2 = (combos(j, 1)==ordinal2);
            connectionVector1 = connectionVectors(combos(j, 1), :, ...
                ordinal1);
            connectionVector2 = connectionVectors(combos(j, 2), :,...
                ordinal2);
            % Compute angles
            angle1 = rad2deg(compute_angle(armVector1, connectionVector1));
            angle2 = rad2deg(compute_angle(armVector2, connectionVector2));
            point_angles(j, :) = [angle1, angle2];
        end

        % Matching function
        matching_value = point_angles(:, 1) + point_angles(:, 2);

        if numel(currentNodeSet)== 2
            % Currently not available
            t=1;
        elseif numel(currentNodeSet)== 3
            [~, argmax] = max(matching_value);
            connectedNodesIndices = combos(argmax, :);
            connectedNodes = currentNodeSet(connectedNodesIndices);
            connectedNodeIDs = G_updated.Nodes.Name(connectedNodes);
            coordPt1 = G_updated.Nodes.Coordinates(connectedNodes(1), :);
            coordPt2 = G_updated.Nodes.Coordinates(connectedNodes(2), :);
            weight = sqrt(sum((coordPt1'-coordPt2').^2));
            G_updated = addedge(G_updated, connectedNodeIDs(1),...
                connectedNodeIDs(2), 1);
        end

    end
end



end
