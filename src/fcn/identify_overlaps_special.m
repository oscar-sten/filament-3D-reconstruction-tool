function G_updated = identify_overlaps_special(current_image,...
    G_input, filter_sigma, structuring_element)


% Note: Disconnected nodes (e.g. through spore occlusion are not considered
G_minimal = removeBridgeNodes(G_input);
degrees = degree(G_minimal);
nodes_degree_3 = find(degrees == 3); % nodes of degree 3

Adj = adjacency(G_minimal);

% Initialize an empty cell array to hold the pairs
node_pairs = [];

pair_edge_weights = [];



% Loop through the nodes of degree 3
for i = 1:length(nodes_degree_3)
    for j = i+1:length(nodes_degree_3)
        % Check if nodes i and j are connected
        if Adj(nodes_degree_3(i),nodes_degree_3(j)) == 1
            % Add the pair to the list
            node_pairs = [node_pairs; nodes_degree_3(i), nodes_degree_3(j)];
        end
    end
end

%% jfslak
n_pairs = size(node_pairs, 1);
margin = 5;

crossings = [];
nodeIDCrossingPairs = [];


for i=1:n_pairs
    % Extract node pair of interest
    node1 = node_pairs(i, 1);
    node2 = node_pairs(i, 2);
    coords1 = G_minimal.Nodes.Coordinates(node1, :);
    coords2 = G_minimal.Nodes.Coordinates(node2, :);
    nodeID1 = G_minimal.Nodes.Name(node1);
    nodeID2 = G_minimal.Nodes.Name(node2);

%     sourceGraph1 = G_minimal.Nodes.SourceGraphIdx(node1);
%     sourceGraph2 = G_minimal.Nodes.SourceGraphIdx(node2);
% 
%     sourceGraphMean = round(mean([sourceGraph1, sourceGraph2]));


%     nodeID1 = G_minimal.Nodes.Name(node1);
%     nodeID2 = G_minimal.Nodes.Name(node2);

    node_coords = [coords1; coords2];

    % Get neighbors
    [P, d] = shortestpath(G_input, nodeID1, nodeID2); % Path between pair

    %     % Check the length of the second shortest path
    %     G_test = rmnode(G_thin, setdiff(P, [nodeID1; nodeID2]));
    %     [~, d2] = shortestpath(G_test, nodeID1, nodeID2); % Path between pair
    %     if ((d/d2)>0.5) && (numel(P)>2)
    %         flag = 1;
    %     else
    %         flag = 0;
    %     end
    flag = testCrossingFeasibility(G_input, nodeID1, nodeID2);
    %     flag = 0;



    nghbrsNode1 = neighbors(G_input, nodeID1);
    nghbrsNode1 = setdiff(nghbrsNode1, P);
    nghbrsNode2 = neighbors(G_input, nodeID2);
    nghbrsNode2 = setdiff(nghbrsNode2, P);

    if numel(nghbrsNode1)==2 && numel(nghbrsNode2)==2
        % Get neighbor coordinates
        nghbr_x_coords = zeros(1, 4);
        nghbr_y_coords = zeros(1, 4);
        for j=1:2
            % Node 1
            nghbrNode = findnode(G_input, nghbrsNode1{j});
            nghbrCoords = G_input.Nodes.Coordinates(nghbrNode, :);
            nghbr_x_coords(2*(j-1)+1) = nghbrCoords(1);
            nghbr_y_coords(2*(j-1)+1) = nghbrCoords(2);
            % Node 2
            nghbrNode = findnode(G_input, nghbrsNode2{j});
            nghbrCoords = G_input.Nodes.Coordinates(nghbrNode, :);
            nghbr_x_coords(2*j) = nghbrCoords(1);
            nghbr_y_coords(2*j) = nghbrCoords(2);
        end

        x_min = min(nghbr_x_coords);
        x_max = max(nghbr_x_coords);
        y_min = min(nghbr_y_coords);
        y_max = max(nghbr_y_coords);

        % Define patch ranges
        x_range = x_min-margin:x_max+margin;
        y_range = y_min-margin:y_max+margin;

        % Remove invalid indices
        [n, m] = size(current_image);
        x_range(x_range<1)=[];
        x_range(x_range>n)=[];
        y_range(y_range<1)=[];
        y_range(y_range>m)=[];

        % Extract and pre-process patch
        patch_test = current_image(x_range, y_range);
        patch_test = mat2gray(patch_test);
%         patch_test = imgaussfilt(patch_test, 0.5);

        % Extract arm vectors
        [armVectors, neighborNodes] = cand_crossing_arm_vectors(G_minimal,...
            G_input, node_pairs(i, :));

        % Compute 2 x 2 cost matrix
        costMatrix = zeros(2, 2);
        for j=1:2
            for k=1:2
                v1 = armVectors(j, :, 1);
                v2 = armVectors(k, :, 2);
                costMatrix(j, k) = 180-rad2deg(compute_angle(v1, v2));
            end
        end

        % Solve optimal matching problem
        % Here it is redundant, since the number of possible solutions is 2.
        matchings = matchpairs(costMatrix, 180);

        pair1 = [armVectors(matchings(:, 1)==1, :, 1)'...
            armVectors(matchings(:, 2)==1, :, 2)'];
        pair2 = [armVectors(matchings(:, 1)==2, :, 1)'...
            armVectors(matchings(:, 2)==2, :, 2)'];

        armVectorPairs = cat(3, pair1, pair2);

        nghbrNodePair1 = [neighborNodes(matchings(:, 1)==1, 1),...
            neighborNodes(matchings(:, 2)==1, 2)];
        nghbrNodePair2 = [neighborNodes(matchings(:, 1)==2, 1),...
            neighborNodes(matchings(:, 2)==2, 2)];

        nghbrNodePairs = cat(3, nghbrNodePair1, nghbrNodePair2);

        if (180-rad2deg(compute_angle(pair1(:, 1), pair1(:, 2)))) > ...
                (180-rad2deg(compute_angle(pair2(:, 1), pair2(:, 2))))
            order = [1, 2];
        else
            order = [2, 1];
        end

        % Assess filaments
        dists = zeros(1, 2);
        dists2 = zeros(1, 2);
        [t, u] = size(patch_test);
        filament_bins = zeros(t, u, 2);
        allNghbrCoords = zeros(2,2,2);
        for j=1:2
            v1 = armVectorPairs(:, 1, order(j));
            v2 = armVectorPairs(:, 2, order(j));
            %             V = v1 - v2;
            %             angle = rad2deg(atan2(V(1), V(2)));
            %             angle = mod(angle+90, 360);
            %             patch_filament_ridge = steerGaussFilterOrder2(patch_test,...
            %                 angle, filter_sigma, 0);
            %             patch_filament_bin = imbinarize(patch_filament_ridge);
            V = v1;
            angle = rad2deg(atan2(V(1), V(2)));
            angle = mod(angle+90, 360);
            patch_filament_ridge = steerGaussFilterOrder2(patch_test,...
                angle, filter_sigma, 0);
            patch_filament_bin1 = imbinarize(patch_filament_ridge);
            patch_filament_bin1 = imdilate(patch_filament_bin1, structuring_element);
            %patch_filament_bin1 = bwmorph(patch_filament_bin1, 'majority');
            V = v2;
            angle = rad2deg(atan2(V(1), V(2)));
            angle = mod(angle+90, 360);
            patch_filament_ridge = steerGaussFilterOrder2(patch_test,...
                angle, filter_sigma, 0);
            patch_filament_bin2 = imbinarize(patch_filament_ridge);
            patch_filament_bin2 = imdilate(patch_filament_bin2, structuring_element);
            %patch_filament_bin2 = bwmorph(patch_filament_bin2, 'majority');
            patch_filament_bin = patch_filament_bin1 & patch_filament_bin2;
            % TEST
            % OBS: If a crossing is very narrow the nodes can result out of
            % the filament, this can be compensated by moving the node
            % sligtly in the direction opposite of the angle perpendicular
            % to the filament

%             patch_filament_bin = bwmorph(patch_filament_bin, 'thicken', 1);
            filament_bins(:, :, j) = patch_filament_bin;

            nbgrNodes = findnode(G_input, string(nghbrNodePairs(:, :,...
                order(j))));
            nghbrCoords = G_input.Nodes.Coordinates(nbgrNodes, :);

            allNghbrCoords(:, :, j) = nghbrCoords;

            in_patch_coords = [nghbrCoords(:, 2)-y_range(1),...
                nghbrCoords(:, 1)-x_range(1)];

            in_patch_node_coords = [node_coords(:, 2)-y_range(1),...
                node_coords(:, 1)-x_range(1)];

            tol = 1;
            n_tol_coords = 2*(tol^2);

            tmp = in_patch_node_coords;
            tmp_idx = sub2ind(size(patch_filament_bin),...
                tmp(:, 2), tmp(:, 1));
            tmp_patch = zeros(size(patch_filament_bin));
            tmp_patch(tmp_idx) = 1;
            tmp_patch = imdilate(tmp_patch, strel('square', tol));
%             [tmp_x, tmp_y] = find(tmp_patch);
%             if numel(tmp_y)==n_tol_coords
%                 in_patch_node_coords = [tmp_y, tmp_x];
%             end
            [lbl, n] = bwlabel(tmp_patch, 8);
            if n==2
                [tmp_x1, tmp_y1] = find(lbl==1);
                [tmp_x2, tmp_y2] = find(lbl==2);
                in_patch_node_coords = [tmp_y1, tmp_x1;
                    tmp_y2, tmp_x2];
%                 % Remove false indices
%                 idx2rm1 = in_patch_node_coords(:, 1)<1;
%                 idx2rm2 = in_patch_node_coords(:, 2)<1;
%                 idx2rm3 = in_patch_node_coords(:, 1)>u-1;
%                 idx2rm4 = in_patch_node_coords(:, 2)>t-1;
%                 in_patch_node_coords((idx2rm1|idx2rm2|idx2rm3|idx2rm4),...
%                     :) = [];

            end
         

            D1 = bwdistgeodesic(patch_filament_bin,...
                in_patch_coords(2,1),...
                in_patch_coords(2, 2));

            dists(j) = D1(in_patch_coords(1, 2), in_patch_coords(1,1));
            if size(in_patch_node_coords, 1)==2
                D2 = bwdistgeodesic(patch_filament_bin,...
                    in_patch_node_coords(2,1),...
                    in_patch_node_coords(2, 2));

                dists2(j) = D2(in_patch_node_coords(1, 2), in_patch_node_coords(1,1));
            else
                
%                 n_tol_coords = size(in_patch_node_coords, 1);

                D2 = bwdistgeodesic(patch_filament_bin,...
                    in_patch_node_coords(...
                    (n_tol_coords/2)+1:n_tol_coords,1),...
                    in_patch_node_coords(...
                   (n_tol_coords/2)+1:n_tol_coords, 2));

                idx = sub2ind(size(D2),...
                    in_patch_node_coords(1:(n_tol_coords/2), 2),...
                    in_patch_node_coords(1:(n_tol_coords/2), 1));

                dists2(j) = min(D2(idx));
            end


        end

%         test = "2855";
%         if nodeID1{1}==test || nodeID2{1}==test
%             figure,
%             imshowpair(filament_bins(:, :, 1), filament_bins(:, :, 2))
%             hold on
%             plot(coords1(2)-y_range(1), coords1(1)-x_range(1), 'b*', 'LineWidth', 2)
%             plot(coords2(2)-y_range(1), coords2(1)-x_range(1), 'b*', 'LineWidth', 2)
%             plot(allNghbrCoords(:, 2, 1)-y_range(1), allNghbrCoords(:, 1, 1)-x_range(1),...
%                 'r*', 'LineWidth',2)
%             plot(allNghbrCoords(:, 2, 2)-y_range(1), allNghbrCoords(:, 1, 2)-x_range(1),...
%                 'y*', 'LineWidth',2)
%             t=1;
%         end

if (((~isnan(dists(1)) && ~isinf(dists(1))) && (~isnan(dists(2)) &&...
        ~isinf(dists(2)))) || ((~isnan(dists2(1)) && ~isinf(dists2(1)))...
        && (~isnan(dists2(2)) &&...
        ~isinf(dists2(2))) &&...
        ((~isnan(dists(1)) && ~isinf(dists(1))) || (~isnan(dists(2)) &&...
        ~isinf(dists(2)))))) && ~flag
    crossings = [crossings; coords1; coords2];
    nodeIDCrossingPairs = [nodeIDCrossingPairs;...
        [str2double(nodeID1{1}), str2double(nodeID2{1})]];

end
    end

end


% x_coords = G_input.Nodes.Coordinates(:, 2);
% y_coords = G_input.Nodes.Coordinates(:, 1);
% 
% figure,
% imshow(current_image)
% hold on
% p=plot(G_input, 'XData', x_coords, 'YData', y_coords, 'LineStyle',...
%     '-','NodeLabel',{},...
%     'EdgeColor', 'r', 'NodeColor', 'b');
% 
% n_test = size(crossings, 1);
% 
% for i=1:2:n_test-1
%     crossings_tmp = crossings(i:i+1, :);
%     plot(crossings_tmp(:, 2), crossings_tmp(:, 1), 'LineWidth', 3)
% end

%% Integration of crossing detections

% TODO: Make sure that the undergoing filament is the one with long
% connection.

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

    % Determine depth of candidates 
    nodePair1 = findnode(G_updated, pair1);
    nodePair2 = findnode(G_updated, pair2);
    depthPair1 = G_updated.Nodes.Z_planes(nodePair1);
    depthPair2 = G_updated.Nodes.Z_planes(nodePair2);

    if mean(depthPair1) > mean(depthPair2)
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
        if numel(mergers) > 2
            mergers = [mergers(1); mergers(end)];
        end
        if numel(mergers) == 2
            nodes2Merge = zeros(1, 2);
            rows2remove = [];
            for j = 1:numel(mergers)
                [x, y] = find(nghbrs2Connect==mergers(j));
                not_y = setdiff([1 2], y);
                corresp = nghbrs2Connect(x, not_y);
                if ~isempty(corresp)
                    nodes2Merge(1, j) = corresp;
                end
                rows2remove = [rows2remove; x];
                %             nghbrs2Connect(x, :) = [];
            end
            if ~any(nodes2Merge==0)
                nghbrs2Connect = [nghbrs2Connect; nodes2Merge];
                nghbrs2Connect(rows2remove, :) = [];
            end
        end
    end
end


%% Merge neighbors

n_pairs2Merge = size(nghbrs2Connect, 1);

for i=1:n_pairs2Merge
    % identify nodes and neighbors
    nodeID1 = nghbrs2Connect(i, 1);
    nodeID2 = nghbrs2Connect(i, 2);
    nghbrs1 = neighbors(G_updated, nodeID1);
    nghbrs2 = neighbors(G_updated, nodeID2);
    % Identify endnodes of edges to remove
    [P, weight] = shortestpath(G_input, nodeID1, nodeID2);
    bp1 = intersect(nghbrs1, P);
    bp2 = intersect(nghbrs2, P);

    % remove edge to BP
    G_updated = rmedge(G_updated, nodeID1, bp1);
    G_updated = rmedge(G_updated, nodeID2, bp2);

    % Add edges to merge
%     [~, weight] = shortestpath(G_thin, nodeID1, nodeID2);
    G_updated = addedge(G_updated, nodeID1, nodeID2, weight);

end




%% Process 4-connected neighbors

% Assuming G is your graph
degrees = degree(G_minimal);
nodes_degree_4 = find(degrees == 4); % nodes of degree 4
nodes_degree_4IDs = G_minimal.Nodes.Name(nodes_degree_4);
perfect_crossing_coords = G_minimal.Nodes.Coordinates(nodes_degree_4, :);

n_4connectedNodes = numel(nodes_degree_4IDs);
for k = 1:n_4connectedNodes

    centerNode = findnode(G_input, nodes_degree_4IDs(k));

    fourNodeNghbrs = neighbors(G_input, centerNode);

    centerCoords = G_input.Nodes.Coordinates(centerNode, :);

    nghbrCoords = G_input.Nodes.Coordinates(fourNodeNghbrs, :);

    armVectors = zeros(4, 2);
    % Here issue arrises at ts 97
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

    pair2mergeIDs = G_input.Nodes.Name(pair2merge);

    centerNodeID = G_input.Nodes.Name(centerNode);

    G_updated = rmedge(G_updated, centerNodeID, pair2mergeIDs(1));
    G_updated = rmedge(G_updated, centerNodeID, pair2mergeIDs(2));

    % Add edges to merge
    [~, weight] = shortestpath(G_input, pair2mergeIDs(1), pair2mergeIDs(2));
    G_updated = addedge(G_updated, pair2mergeIDs(1), pair2mergeIDs(2), weight);
end
