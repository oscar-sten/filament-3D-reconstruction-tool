function [nodeIDCrossingPairs, crossings] = identify_overlaps(...
    current_image, G_input, G_minimal, filter_sigma, ...
    side_margin, overlap_displacement_tolerance)


if overlap_displacement_tolerance>0
    tol = 2*(overlap_displacement_tolerance)+1;
end


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
margin = side_margin;

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
    node_coords = [coords1; coords2];

    % Get neighbors
    [P, ~] = shortestpath(G_input, nodeID1, nodeID2); % Path between pair
    
    % Simple sanity check. In cases where there are loop-like structures.
    flag = testCrossingFeasibility(G_input, nodeID1, nodeID2);

    nghbrsNode1 = neighbors(G_input, nodeID1);
    nghbrsNode1 = setdiff(nghbrsNode1, P);
    nghbrsNode2 = neighbors(G_input, nodeID2);
    nghbrsNode2 = setdiff(nghbrsNode2, P);

    nghbrsNodesBoth = unique([nghbrsNode1; nghbrsNode2]);

    if (numel(nghbrsNode1)==2 && numel(nghbrsNode2)==2) &&...
            (numel(nghbrsNodesBoth)==4)
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
        
        % Make sure that no important coordinate is outside the patch
        x_min = min([nghbr_x_coords, coords1(1), coords2(1)]);
        x_max = max([nghbr_x_coords, coords1(1), coords2(1)]);
        y_min = min([nghbr_y_coords, coords1(2), coords2(2)]);
        y_max = max([nghbr_y_coords, coords1(2), coords2(2)]);

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
            V = v1;
            angle = rad2deg(atan2(V(1), V(2)));
            angle = mod(angle+90, 360);
            patch_filament_ridge = steerGaussFilterOrder2(patch_test,...
                angle, filter_sigma, 0);
            patch_filament_bin1 = imbinarize(patch_filament_ridge);
            V = v2;
            angle = rad2deg(atan2(V(1), V(2)));
            angle = mod(angle+90, 360);
            patch_filament_ridge = steerGaussFilterOrder2(patch_test,...
                angle, filter_sigma, 0);
            patch_filament_bin2 = imbinarize(patch_filament_ridge);
            patch_filament_bin = patch_filament_bin1 & patch_filament_bin2;

            filament_bins(:, :, j) = patch_filament_bin;

            nbgrNodes = findnode(G_input, string(nghbrNodePairs(:, :,...
                order(j))));
            nghbrCoords = G_input.Nodes.Coordinates(nbgrNodes, :);

            allNghbrCoords(:, :, j) = nghbrCoords;
            
            % Current
            in_patch_coords = [nghbrCoords(:, 2)-y_range(1),...
                nghbrCoords(:, 1)-x_range(1)];

            in_patch_node_coords = [node_coords(:, 2)-y_range(1),...
                node_coords(:, 1)-x_range(1)];

%             % Test
%             in_patch_coords = [nghbrCoords(:, 2)-y_range(1)+1,...
%                 nghbrCoords(:, 1)-x_range(1)+1];
%             in_patch_node_coords = [node_coords(:, 2)-y_range(1)+1,...
%                 node_coords(:, 1)-x_range(1)+1];
            
            if overlap_displacement_tolerance>0
                
                n_tol_coords = 2*(tol^2);

                tmp = in_patch_node_coords;
                tmp_idx = sub2ind(size(patch_filament_bin),...
                    tmp(:, 2), tmp(:, 1));
                tmp_patch = zeros(size(patch_filament_bin));
                tmp_patch(tmp_idx) = 1;
                tmp_patch = imdilate(tmp_patch, strel('square', tol));

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
            end
         

            D1 = bwdistgeodesic(patch_filament_bin,...
                in_patch_coords(2,1),...
                in_patch_coords(2, 2));


            dists(j) = D1(in_patch_coords(1, 2), in_patch_coords(1,1));
            if size(in_patch_node_coords, 1)==2
                D2 = bwdistgeodesic(patch_filament_bin,...
                    in_patch_node_coords(2,1),...
                    in_patch_node_coords(2, 2));
                try
                dists2(j) = D2(in_patch_node_coords(1, 2), in_patch_node_coords(1,1));
                catch
                    figure,
                    imshowpair(filament_bins(:, :, 1), filament_bins(:, :, 2))
                    hold on
                    plot(in_patch_node_coords(2, :), in_patch_node_coords(1, :), '*')
                    t=1;
                end
            else

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
        
%         figure,
%         imshowpair(filament_bins(:, :, 1), filament_bins(:, :, 2))
%         hold on
%         plot(in_patch_node_coords(2, :), in_patch_node_coords(1, :), '*')
%         t=1;

        if ((~isnan(dists(1)) && ~isinf(dists(1))) && (~isnan(dists(2)) &&...
                ~isinf(dists(2)))) || ((~isnan(dists2(1)) && ~isinf(dists2(1)))...
                && (~isnan(dists2(2)) &&...
                ~isinf(dists2(2))) &&...
                ((~isnan(dists(1)) && ~isinf(dists(1))) || (~isnan(dists(2)) &&...
                ~isinf(dists(2))))) && ~flag
            crossings = [crossings; coords1; coords2];
            nodeIDCrossingPairs = [nodeIDCrossingPairs;...
                [str2double(nodeID1{1}), str2double(nodeID2{1})]];

        end
    end

end



end