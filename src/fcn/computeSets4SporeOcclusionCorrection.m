function [nodeSets, armVectorSet, connectionVectorSet] =...
    computeSets4SporeOcclusionCorrection(G_updated, idSporeMask,...
    allSporeMasks)


% Requirements: idSporeMask, all 1-deg nodes.

deg1Nodes = find(degree(G_updated)==1);

deg1Coords = G_updated.Nodes.Coordinates(deg1Nodes, :);

ind = sub2ind(size(idSporeMask), deg1Coords(:, 1), deg1Coords(:, 2));

% Create expanded spore mask
n_spore_masks = size(allSporeMasks, 3);
idSporeMask = zeros(size(idSporeMask));

for i=1:n_spore_masks
    idSporeMask = idSporeMask | allSporeMasks(:, :, i);
end

idSporeMask = bwlabel(idSporeMask);

nodeSets = {};

for i = 1:max(idSporeMask(:))
    mask = (idSporeMask == i);
    nodeSet = deg1Nodes(mask(ind));
    if ~isempty(nodeSet)
        nodeSets{end+1} = nodeSet;
    end

end


armVectorSet = {};
connectionVectorSet = {};
for i=1:numel(nodeSets)
    currentNodeSet = nodeSets{i};
    if numel(currentNodeSet) > 1
        armVectors = zeros(numel(currentNodeSet), 2);

        for j=1:numel(currentNodeSet)
            %             connectionNode = G_updated.Nodes.Name(currentNodeSet(j));
            connectionNode = currentNodeSet(j);
            nghbr = neighbors(G_updated, connectionNode);
            connectionCoords = G_updated.Nodes.Coordinates(...
                connectionNode, :);
            nghbrCoords = G_updated.Nodes.Coordinates(...
                nghbr, :);
          
            degNghbr = degree(G_updated, nghbr);
            if degNghbr == 2
                secondaryNeighbors = neighbors(G_updated, nghbr);
                secondaryNeighbor = setdiff(secondaryNeighbors, ...
                    connectionNode);
                sec_nghbr = secondaryNeighbor;
                nghbrCoords2 = G_updated.Nodes.Coordinates(...
                    sec_nghbr, :);
                armVector = nghbrCoords2 - nghbrCoords;
            else
                armVector = nghbrCoords - connectionCoords;
            end

            armVectors(j, :) = armVector;
        end

        connectionVectors = zeros(numel(currentNodeSet), 2,...
            numel(currentNodeSet)-1);
        for j=1:numel(currentNodeSet)
            currentNode = currentNodeSet(j);
            otherNodes = setdiff(currentNodeSet, currentNode);
            currentNodeCoords = G_updated.Nodes.Coordinates(...
                currentNode, :);
            for k=1:numel(currentNodeSet)-1
                currentOtherNodeCoords = G_updated.Nodes.Coordinates(...
                    otherNodes(k), :);
                connectionVector = currentOtherNodeCoords...
                    - currentNodeCoords;
                connectionVectors(j, :, k) = connectionVector;
            end
        end
    else
        armVectors = [];
        connectionVectors = [];
    end
    armVectorSet{end+1} = armVectors;
    connectionVectorSet{end+1} = connectionVectors;
end


end