function G_merged = mergeGraphs(G_thin_list, imageMetaData,....
    depthParams, pruneAngleThreshold)

% Image stack meta data
data_path_prefix = imageMetaData.dataPathPrefix;
folder = imageMetaData.folder;
files = imageMetaData.files;
selectedPlane = imageMetaData.selectedPlane;
planeNumber = imageMetaData.planeNumber;


% Depth estimation parameters
patch_range = depthParams.patchRange;
FM_type = depthParams.FMType;
costUnmatched = depthParams.costUnmatched;


% Merge graphs iteratively
G_best_so_far = [];
debug_mode = 0;

for i = 1:numel(G_thin_list)

    if isempty(G_best_so_far)
        % If no previous graph has been processed
        G_best_so_far = G_thin_list{i};
        G_best_so_far.Nodes.OriginalDegree = degree(G_best_so_far);
        G_best_so_far.Nodes.SourceGraphIdx = NaN(size(G_best_so_far.Nodes, 1), 1);

        % Compute
        n_threads = size(G_best_so_far.Nodes, 1);
        current_point_threads = G_best_so_far.Nodes.Coordinates;
        focusValues = zeros(n_threads, 1);
        img_path = strcat(data_path_prefix, folder,...
            files(planeNumber(i)).name);
        tmp_image_color = imread(img_path);

        if size(tmp_image_color, 3)==2
            tmp_image = rgb2gray(tmp_image_color);
        else
            tmp_image = tmp_image_color;
        end

        for ii=1:n_threads
            coords = current_point_threads(ii, :);
            [patch, ~, ~] =...
                extractPatch(tmp_image, coords, patch_range);
            FM = fmeasure(patch, FM_type);
            focusValues(ii) = FM;
        end

        G_best_so_far.Nodes.Node_FM_value = focusValues;
        G_best_so_far.Nodes.FM_list = cell(size(G_best_so_far.Nodes, 1), 1);
        G_best_so_far.Nodes.FM_list_planes = ...
            cell(size(G_best_so_far.Nodes, 1), 1);
        for ii=1:n_threads
            G_best_so_far.Nodes.FM_list{ii} = focusValues(ii);
            G_best_so_far.Nodes.FM_list_planes{ii} = planeNumber(i);
        end
    else
        % If a previous plane has been processed
        % Formulate the matching problem between nodes in the 
        % 'best-so-far'-graph and the "new" graph.
        graph1 = G_best_so_far;
        graph1_coords = graph1.Nodes.Coordinates;
        graph1_IDs = graph1.Nodes.Name;

        graph2 = G_thin_list{i};

        graph2_coords = graph2.Nodes.Coordinates;

        graph2_IDs = graph2.Nodes.Name;

        n_nodes_graph1 = size(graph1_coords, 1);

        n_nodes_graph2 = size(graph2_coords, 1);

        % Define cost matrix
        costMatrix = zeros(n_nodes_graph1, n_nodes_graph2);

        for ii = 1:n_nodes_graph1
            for iii = 1:n_nodes_graph2
                costMatrix(ii, iii) = norm(graph1_coords(ii, :)-...
                    graph2_coords(iii, :));
            end
        end

        % Solve matching problem
        [matchings, unassignedRows, unassignedCols] = matchpairs(costMatrix,...
            costUnmatched);

        first_matched = graph1_coords(matchings(:, 1), :);
        second_matched = graph2_coords(matchings(:, 2), :);

        first_unmatched = graph1_coords(unassignedRows,:);
        second_unmatched = graph2_coords(unassignedCols, :);

        % Log Node IDs
        first_matched_IDs = cellfun(@string, graph1_IDs(matchings(:, 1)));
        second_matched_IDs = cellfun(@string, graph2_IDs(matchings(:, 2)));

        first_unmatched_IDs = cellfun(@string, graph1_IDs(unassignedRows,:));
        second_unmatched_IDs = cellfun(@string, graph2_IDs(unassignedCols, :));

        nodeID_threads = [first_matched_IDs, second_matched_IDs];

        current_point_threads = [first_matched, second_matched];


        single_plane_nodes = [[ones(size(...
            first_unmatched_IDs, 1), 1)*(1),...
            first_unmatched_IDs]; [ones(size(...
            second_unmatched_IDs, 1), 1)*(2),...
            second_unmatched_IDs]];
        single_plane_points_list = [[(1)*ones(size(...
            first_unmatched, 1), 1) first_unmatched];...
            [(2)*ones(size(...
            second_unmatched, 1), 1) second_unmatched]];


        % Make comparison for all point threads

        [n_threads, ~] = size(nodeID_threads);

        % Comupute focus values for new graph
        focusValues = NaN(n_threads, 1);

        img_path = strcat(data_path_prefix, folder,...
            files(planeNumber(i)).name);
        tmp_image_color = imread(img_path);
        
        if size(tmp_image_color, 3)==2
            tmp_image = rgb2gray(tmp_image_color);
        else
            tmp_image = tmp_image_color;
        end
        
        for ii=1:n_threads
            coords = current_point_threads(ii, 3:4);
            [patch, ~, ~] =...
                extractPatch(tmp_image, coords, patch_range);
            FM = fmeasure(patch, FM_type);
            focusValues(ii) = FM;
        end

        idx_old_matched = findnode(G_best_so_far, nodeID_threads(:, 1));
        FM_values_best_so_far = G_best_so_far.Nodes.Node_FM_value(...
            idx_old_matched);


        % Compute single plane node focus values
        single_node_FM_values = [];
        new_single_points_idx = find(single_plane_points_list(:, 1)==2);
        if ~isempty(new_single_points_idx)
            new_single_points = single_plane_points_list(...
                new_single_points_idx, 2:3);
%             img_path = strcat(data_path_prefix, folder,...
%                 files(selectedPlane(i)).name);
            img_path = strcat(data_path_prefix, folder,...
                files(i).name);

            tmp_image_color = imread(img_path);
            tmp_image = rgb2gray(tmp_image_color);

            for ii=1:size(new_single_points, 1)
                current_coords = new_single_points(ii, :);
                [patch, ~, ~] =...
                    extractPatch(tmp_image, current_coords, patch_range);
                FM = fmeasure(patch, FM_type);
                single_node_FM_values = [single_node_FM_values; FM];
            end

        end




        % Compare
        focusValues_cmp = [FM_values_best_so_far, focusValues];

        [FMmax, planeIndexWithMaxFM] = max(focusValues_cmp , [], 2);

% 
%         planeWithMaxFM = selectedPlane(planeIndexWithMaxFM)';

        bestCoordinates = zeros(size(planeIndexWithMaxFM, 1), 2);
        bestNodes = zeros(size(planeIndexWithMaxFM, 1), 1);

        for ii=1:n_threads
            bestCoordinates(ii, :) = current_point_threads(ii, ...
                2*planeIndexWithMaxFM(ii)-1:2*planeIndexWithMaxFM(ii));
            bestNodes(ii) = nodeID_threads(ii, planeIndexWithMaxFM(ii));
        end


        % Reconstruct optimal graph
        n_new_nodes = size(planeIndexWithMaxFM, 1)+size(...
            single_plane_points_list, 1);

        new_node_list = string((1:n_new_nodes)');


        G_test = graph(zeros(n_new_nodes, n_new_nodes), new_node_list);

        if isempty(single_plane_points_list)
            bestCoordinatesAll = bestCoordinates;
            bestNodesAll = bestNodes;
            sourceThread = planeIndexWithMaxFM;
        else
            bestCoordinatesAll = [bestCoordinates; single_plane_points_list(:, 2:3)];
            bestNodesAll = [bestNodes; single_plane_nodes(:, 2)];
            sourceThread = [planeIndexWithMaxFM; single_plane_points_list(:, 1)];
        end


        G_test.Nodes.SourceGraphIdx = NaN(size(G_test.Nodes, 1), 1);
        G_test.Nodes.Node_FM_value = NaN(size(G_test.Nodes, 1), 1);

        G_test.Nodes.FM_list = cell(size(G_test.Nodes, 1), 1);

        node_correspondance = zeros(size(bestCoordinatesAll));

        for ii=1:n_new_nodes
            currentNode = new_node_list(ii);
            best_thread = sourceThread(ii);
            if best_thread < 2
                current_graph = G_best_so_far;
            else
                current_graph = G_thin_list{i};
            end
            corresponding_node = bestNodesAll(ii);

            % Store the computed FM values and corresponding plane indices
            % for each node.
            % First step: identify type: matched, former (in previous but
            % not present) or new.
            % If matched: retive old info and add new.
            if ii<=n_threads
                % Identify previous node ID
                previous_idx = matchings(ii, 1);
                % Log FM value
                previous_FM_list = G_best_so_far.Nodes.FM_list{...
                    previous_idx};
                new_FM_value = focusValues(ii);
                G_test.Nodes.FM_list{ii} = [previous_FM_list new_FM_value];
                % Log plane number
                previous_FM_list_planes = ...
                    G_best_so_far.Nodes.FM_list_planes{...
                    previous_idx};
                G_test.Nodes.FM_list_planes{ii} = [previous_FM_list_planes...
                    planeNumber(i)];
            % If former: simpy add all old info.
            elseif single_plane_points_list(ii-n_threads, 1) == 1
                previous_idx = unassignedRows(ii-n_threads, 1);
                previous_FM_list = G_best_so_far.Nodes.FM_list{...
                    previous_idx};
                G_test.Nodes.FM_list{ii} = previous_FM_list;
                previous_FM_list_planes =...
                    G_best_so_far.Nodes.FM_list_planes{...
                    previous_idx};
                G_test.Nodes.FM_list_planes{ii} = previous_FM_list_planes;
            % If new: simply add new values.
            elseif single_plane_points_list(ii-n_threads, 1) == 2
                G_test.Nodes.FM_list{ii} = single_node_FM_values(...
                    ii-(n_threads+numel(unassignedRows)));
                G_test.Nodes.FM_list_planes{ii} = planeNumber(i);
            else
                error("Invalid plane index")
            end



            if ~isnan(G_best_so_far.Nodes.SourceGraphIdx)
                if best_thread < 2
                    idx = findnode(current_graph, corresponding_node);
                    G_test.Nodes.SourceGraphIdx(ii) = ...
                        current_graph.Nodes.SourceGraphIdx(idx);
                    G_test.Nodes.Node_FM_value(ii) = ...
                        current_graph.Nodes.Node_FM_value(idx);

                elseif ii <= numel(bestNodes)
                    G_test.Nodes.SourceGraphIdx(ii) = i;
                    G_test.Nodes.Node_FM_value(ii) = FMmax(ii);

                else
                    G_test.Nodes.SourceGraphIdx(ii) = i;
                    G_test.Nodes.Node_FM_value(ii) =...
                        single_node_FM_values(ii-numel(bestNodes)-...
                        numel(find(single_plane_points_list(:, 1)==1)));

                end
            end

            [~,nid] = outedges(current_graph, string(corresponding_node));
            node_correspondance(ii, :) = [currentNode, string(corresponding_node)];

            for iii=1:numel(nid)
                if ii <= numel(bestNodes)
                    % If the node has been matched
                    correspIdx = find(string(nodeID_threads(:,...
                        sourceThread(ii)))==nid(iii));
                    idxOut = findedge(current_graph,...
                        string(corresponding_node), nid(iii));
                    weight = current_graph.Edges.Weight(idxOut);

                else
                    % If the node is present only in a single plane
                    correspIdx = find(string(best_thread)==single_plane_nodes(:, 1) & ...
                        single_plane_nodes(:, 2)==nid(iii));

                    if sum(string(best_thread)==single_plane_nodes(:, 1) & ...
                            single_plane_nodes(:, 2)==nid(iii)) == 0
                        % Check is connection in current threads exist
                        test = find(string(nodeID_threads(:, sourceThread(ii)))==nid(iii));
%                         correspNode = nodeID_threads(test, sourceThread(test));
                        correspIdx = test;
                    else
                        correspIdx = correspIdx + numel(bestNodes);
                    end

                    idxOut = findedge(current_graph, string(corresponding_node), nid(iii));
                    weight = current_graph.Edges.Weight(idxOut);
                end

                % Add edge
                G_test = addedge(G_test, currentNode,...
                    G_test.Nodes.Name(correspIdx), weight);

            end

        end

        if isnan(G_best_so_far.Nodes.SourceGraphIdx)
            sourceGraph = zeros(size(G_test.Nodes, 1), 1);
            bestFMs = zeros(size(G_test.Nodes, 1), 1);
            sourceGraph(1:numel(bestNodes)) = planeIndexWithMaxFM;
            if ~isempty(single_plane_nodes)
                sourceGraph(numel(bestNodes)+find(...
                    single_plane_nodes(:, 1)=="1")) = i-1;
                sourceGraph(numel(bestNodes)+find(...
                    single_plane_nodes(:, 1)=="2")) = i;
            end
            % Update best FMs
            bestFMs(1:numel(bestNodes)) = FMmax;
            if ~isempty(single_plane_nodes)
                bestFMs(numel(bestNodes)+find(...
                    single_plane_nodes(:, 1)=="1")) = ...
                    G_best_so_far.Nodes.Node_FM_value(numel(bestNodes)+find(...
                    single_plane_nodes(:, 1)=="1"));
                bestFMs(numel(bestNodes)+find(...
                    single_plane_nodes(:, 1)=="2")) = single_node_FM_values;
            end




            G_test.Nodes.SourceGraphIdx = sourceGraph;

            G_test.Nodes.Node_FM_value = bestFMs;
        end
        sourceGraph = G_test.Nodes.SourceGraphIdx;

        optimal_z_planes = planeNumber(sourceGraph);

        G_test.Nodes.Z_planes = optimal_z_planes';


        G_test = simplify(G_test);

        G_test.Nodes.Coordinates = [bestCoordinatesAll(:, 1),...
            bestCoordinatesAll(:, 2)];


        % All cycles instead of cyclebasis is correct, cycle basis produces errors.
        [cycles,edgecycles] = allcycles(G_test);
        nodeDegsMerged = degree(G_test);
%         nodeIDsMerged = G_test.Nodes.Name;

        % Identify the nodes that have a different (higher) degree with respect to
        % the original node in the origin graph.

        % Make lists of node IDs in new graph and in origin graph
        original_nodes_degrees = zeros(numel(nodeDegsMerged), 1);
        diffDegNodes = [];

        for ii=1:numel(nodeDegsMerged)
            currentNode = G_test.Nodes.Name(ii);
            currentNodeDeg = degree(G_test, currentNode);
            originNode = node_correspondance(ii, 2);
            if sourceThread(ii)==1
                idx = findnode(G_best_so_far, string(originNode));
                originNodeDeg = G_best_so_far.Nodes.OriginalDegree(idx);
            else
                origin_graph = G_thin_list{i};
                originNodeDeg = degree(origin_graph, string(originNode));
            end

            original_nodes_degrees(ii) = originNodeDeg;
            if currentNodeDeg ~= originNodeDeg
                diffDegNodes = [diffDegNodes; currentNode];
            end
        end


        G_test.Nodes.OriginalDegree = original_nodes_degrees;

        % Identify separated paths that should be together

        % START WITH COMBINATIONS BY CYCLE!
        setsOfPaths = {};
        setsOfEdgePaths = {};
        endNodes = [];

        cand_set = {};
        corresp_cycle = [];
        in_a_cycle = [];

        for ii=1:numel(cycles)
            common = intersect(cycles{ii}, diffDegNodes);
            if numel(common)>1
                cand_set{end+1} = common;
                corresp_cycle = [corresp_cycle; ii];
            end
            in_a_cycle = [in_a_cycle; common];
        end

        for ii=1:numel(corresp_cycle)
            currentCycle = cycles{corresp_cycle(ii)};
            currentEdgeCycle = edgecycles{corresp_cycle(ii)};
            currentSubGraph = subgraph(G_test, currentCycle);
            graphTmp = G_test;
            graphTmp = rmedge(graphTmp, setdiff(1:numel(graphTmp.Edges.Weight),...
                currentEdgeCycle));
            currentSet = cand_set{ii};
            diffDegNodes_test = cellfun(@str2double, currentSet);
            diffDeg_combos = nchoosek(diffDegNodes_test, 2);
            n_combos = size(diffDeg_combos, 1);
            path_node_set = [];
            dist_combos = zeros(n_combos, 2);
            for iii=1:n_combos
                % In a cycle in a graph, between each pair of nodes, there are two
                % paths.
                [~,~,edgepath] = shortestpath(graphTmp,...
                    string(diffDeg_combos(iii, 1)), string(diffDeg_combos(iii, 2)));
                otherEdgePath = setdiff(currentEdgeCycle, edgepath);
                dist_combos(iii, 1) = sum(graphTmp.Edges.Weight(edgepath));
                graphTmp2 = rmedge(graphTmp, edgepath);
                dist_combos(iii, 2) = sum(graphTmp2.Edges.Weight);

            end

            quote = dist_combos(:, 2)./dist_combos(:, 1);
            % The most likely pair is the one with most similar distances
            [~, idx] = min(quote);
            candidate_nodes = diffDeg_combos(idx, :);
            angles = zeros(numel(candidate_nodes), 1);
            for iii=1:numel(candidate_nodes)
                % Compute direction vectors for edges in edge cycle.
                [eid,nid] = outedges(graphTmp, string(candidate_nodes(iii)));
                %         node_edges_in_cycle = intersect(eid, currentEdgeCycle);
                current_endNodes = graphTmp.Edges.EndNodes(eid, :);
                source_node = candidate_nodes(iii);
                target_nodes = setdiff(cellfun(@string, current_endNodes),...
                    string(source_node));
                idx_src = findnode(graphTmp, source_node);
                source_coords = graphTmp.Nodes.Coordinates(idx_src, :);
                idx_target = findnode(graphTmp, target_nodes);
                target_coords = graphTmp.Nodes.Coordinates(idx_target, :);
                direction_vectors = source_coords-target_coords;
                phi = compute_angle(direction_vectors(1, :)',...
                    direction_vectors(2, :)');
                phi_deg = rad2deg(phi);
                angles(iii) = phi_deg;
            end


            if numel(currentCycle)<=3
                [nodePath, ~, edgepath] = shortestpath(graphTmp,...
                    string(candidate_nodes(1)),...
                    string(candidate_nodes(2)));

                graphTmp2 = rmedge(graphTmp, edgepath);

                [otherNodePath, ~, otherEdgePath] = shortestpath(graphTmp2,...
                    string(candidate_nodes(1)),...
                    string(candidate_nodes(2)));

                setsOfPaths{end+1} = {nodePath; otherNodePath};
                setsOfEdgePaths{end+1} = {edgepath; otherEdgePath};
                endNodes = [endNodes;...
                    string(candidate_nodes)];

            elseif all(angles<pruneAngleThreshold)
                [nodePath, ~, edgepath] = shortestpath(graphTmp,...
                    string(candidate_nodes(1)),...
                    string(candidate_nodes(2)));

                graphTmp2 = rmedge(graphTmp, edgepath);

                [otherNodePath, ~, otherEdgePath] = shortestpath(graphTmp2,...
                    string(candidate_nodes(1)),...
                    string(candidate_nodes(2)));


                setsOfPaths{end+1} = {nodePath; otherNodePath};
                setsOfEdgePaths{end+1} = {edgepath; otherEdgePath};
                endNodes = [endNodes;...
                    string(candidate_nodes)];
            end

        end


        % Identify tails using direction strategy
        % In this step we might investigate all non-matched endnodes, on the
        % off-chance that the start of the tails is on a cycle.

        % A tail may beginn inside a cycle!
        idx = findnode(G_test, diffDegNodes);
        degreeDiff = degree(G_test, diffDegNodes)-G_test.Nodes.OriginalDegree(idx);
        if ~isempty(endNodes)
            not_in_a_cycle = setdiff(diffDegNodes, reshape(endNodes,...
                [numel(endNodes), 1]));
        else
            not_in_a_cycle = [];
        end

        deg_not_in_cycle = degree(G_test, not_in_a_cycle);
        not_in_a_cycle(deg_not_in_cycle<3) = [];



        % Step 1: Identify the edges with smallest angle between them.
        % Step 2: Check that the angle is resonably small.
        % Step 3: Check that the edges are not in any cycle.
        not_in_a_cycle = [];

        for ii=1:numel(diffDegNodes)
            current_node = diffDegNodes(ii);
            [eid,nidStart] = outedges(G_test, current_node);

            idx_base = findnode(G_test, current_node);
            base_coords = G_test.Nodes.Coordinates(idx_base, :);
            cand_tail_V = zeros(size(eid, 1), 2);
            if numel(nidStart) > 1
                for iii=1:numel(nidStart)
                    idx_goal = findnode(G_test, nidStart(iii));
                    goal_coords = G_test.Nodes.Coordinates(idx_goal, :);
                    direction_vector = base_coords-goal_coords;
                    cand_tail_V(iii, :) = direction_vector;
                end
                % Compute angles between candidate vectors
                combos = nchoosek(1:numel(nidStart), 2);
                cand_tail_phi = zeros(size(combos, 1), 1);

                for iii=1:size(combos, 1)
                    idx = combos(iii, :);
                    direction_vectors = cand_tail_V(idx', :);
                    phi = compute_angle(direction_vectors(1, :)',...
                        direction_vectors(2, :)');
                    phi_deg = rad2deg(phi);
                    cand_tail_phi(iii) = phi_deg;
                end

                % In most cases, the smallest angle is very small
                % Determine if two tails are close to each other
                idx = find(cand_tail_phi < pruneAngleThreshold);
                cand = eid(combos(idx, :));

                for iii=1:numel(idx)
                    common_all = [];
                    for iv=1:numel(edgecycles)
                        if numel(idx)>1
                            common = intersect(edgecycles{iv}, cand(iii, :));
                        else
                            common = intersect(edgecycles{iv}, cand);
                        end
                        if ~isempty(common)
                            common_all = [common_all; ...
                                reshape(common, [numel(common), 1])];
                        end

                    end

                    if numel(unique(common_all))<2
                        not_in_a_cycle = [not_in_a_cycle; current_node];
                    end
                end
            end
        end

        identifiedTails = {};

        if ~isempty(not_in_a_cycle)

            for ii=1:numel(not_in_a_cycle)
                current_remainingDiffDegNode = not_in_a_cycle(ii);
                [eid,nidStart] = outedges(G_test, current_remainingDiffDegNode);

                % Investigate if there are two nodes with small angle between them
                idx_base = findnode(G_test, current_remainingDiffDegNode);
                base_coords = G_test.Nodes.Coordinates(idx_base, :);
                cand_tail_V = zeros(size(nidStart, 1), 2);

                for iii=1:numel(nidStart)
                    idx_goal = findnode(G_test, nidStart(iii));
                    goal_coords = G_test.Nodes.Coordinates(idx_goal, :);
                    direction_vector = base_coords-goal_coords;
                    cand_tail_V(iii, :) = direction_vector;
                end

                % Compute angles between candidate vectors
                combos = nchoosek(1:numel(nidStart), 2);
                cand_tail_phi = zeros(size(combos, 1), 1);

                for iii=1:size(combos, 1)
                    idx = combos(iii, :);
                    direction_vectors = cand_tail_V(idx', :);
                    phi = compute_angle(direction_vectors(1, :)',...
                        direction_vectors(2, :)');
                    phi_deg = rad2deg(phi);
                    cand_tail_phi(iii) = phi_deg;
                end

                % Determine if two tails are close to each other
                idx = (cand_tail_phi < pruneAngleThreshold);
                tailStarts = nidStart(combos(idx, :));

                pathsAway = {};


                for iii=1:numel(tailStarts)
                    current_tail = [current_remainingDiffDegNode; tailStarts(iii)];
                    current_Node = current_tail(end);
                    [eid,nid] = outedges(G_test, current_Node);
                    diff = setdiff(nid, current_tail);
                    if numel(diff)>1
                        diff=[];
                        current_tail = [];
                    end
                    while ~isempty(diff) && ~isempty(current_Node)
                        current_tail = [current_tail; diff];
                        current_Node = current_tail(end);
                        [eid,nid] = outedges(G_test, current_Node);
                        diff = setdiff(nid, current_tail);
                        if numel(diff)>1
                            % If we arrive to a cross-road.
                            % 1. Check original degree
                            idx = findnode(G_test, current_Node);
                            originalDeg = G_test.Nodes.OriginalDegree(idx);
                            diff=[];
                            current_tail = [];

                        end
                    end
                    if ~isempty(current_tail)
                        tail_tmp = cellfun(@string, current_tail);
                        %             tail_tmp = setdiff(tail_tmp, current_remainingDiffDegNode);
                        pathsAway{end+1} = tail_tmp;
                    end
                end
                % Add tails
                identifiedTails{end+1} = pathsAway;
            end


        end



        if debug_mode
            nodeColors = ['r', 'm'];
            img_path = strcat(data_path_prefix, folder, files(15).name);
            Img = imread(img_path);
            figure,
            imshow(Img)
            hold on
            p=plot(G_test, 'XData', bestCoordinatesAll(:, 2),...
                'YData',bestCoordinatesAll(:, 1),...
                'LineStyle',...
                '-','NodeLabel',{},...
                'EdgeColor', 'w', 'NodeColor', 'r', 'MarkerSize', 3, 'LineWidth', 2);

            for ii=1:numel(setsOfPaths)
                currentPathPair = setsOfPaths{ii};
                outerNodes = endNodes(ii, :);
                for iii=1:numel(currentPathPair)
                    currentNodes = currentPathPair{iii};
                    %         currentNodes = currentNodes{1};
                    pathNodes = setdiff(currentNodes, outerNodes);
                    highlight(p, pathNodes, 'NodeColor', nodeColors(iii))

                end
                highlight(p, outerNodes, 'NodeColor', 'b')

            end

            highlight(p, not_in_a_cycle, 'NodeColor', 'cyan')

            nodeColors = ['g', 'y'];


            for ii=1:numel(identifiedTails)
                current_tail = identifiedTails{ii};
                current_remainingDiffDegNode = not_in_a_cycle(ii);
                for iii=1:numel(current_tail)
                    current_sub_tail = current_tail{iii};
                    current_sub_tail = setdiff(current_sub_tail,...
                        current_remainingDiffDegNode);
                    highlight(p, current_sub_tail, 'NodeColor', nodeColors(iii));
                end
            end
            title('Identified cyles and tails')
        end

        G_merged = G_test;

        node_degs_test = degree(G_test);

        % Remove faulty paths in cycles
        for ii=1:numel(setsOfPaths)
            paths = setsOfPaths{ii};
            edgePaths = setsOfEdgePaths{ii};
            pathLengths = cellfun(@numel, paths);
            currentPair = endNodes(ii, :);
            if any(pathLengths==2)

                idx = (pathLengths==2);

                if sum(idx)==1
                    longer_path = paths{~idx};
                    shorter_path = paths{idx};

                    test = degree(G_test, setdiff(longer_path, ...
                        shorter_path))>2;
                    if sum(test)==1
                        toRemove = shorter_path;
                    else
                        toRemove = longer_path;
                    end

                else
                    toRemove = paths{idx};
                end

            else
                path_focus_values = zeros(numel(paths), 1);
                branch_presence = zeros(numel(paths), 1);
                for iii=1:numel(paths)
                    otherNodes = setdiff(paths{iii}, currentPair');
                    idx = findnode(G_test, otherNodes);
                    path_focus_values(iii) = median(G_test.Nodes.Node_FM_value(idx));
                    path_node_degs = node_degs_test(idx);
                    deg3 = (path_node_degs>2);
                    % Check if there is a branch
                    if  ~isempty(G_test.Nodes.Name(idx(deg3)))
                        if isempty(intersect(G_test.Nodes.Name(idx(deg3)),...
                                diffDegNodes)) || ...
                                G_test.Nodes.OriginalDegree(idx(deg3))>2
                            branch_presence(iii) = 1;
                        end
                    end

                end

                % All paths with branch are kept (i.e., if both are with
                % branch neither is removed). If both are without branch,
                % selection is based on FM value.
                if sum(branch_presence)==1
                    idx_no_branch = find(~branch_presence);
                    toRemove = [];
                    for iii=1:numel(idx_no_branch)
                        toRemove = [toRemove, paths{idx_no_branch(iii)}];
                    end
                elseif sum(branch_presence)==0
                    [~, idx] = min(path_focus_values);
                    toRemove = paths{idx};
                elseif sum(branch_presence)==2
                    % This is a special case that occurs when a crossover
                    % point in the previous plane results as a pair of
                    % branch points, and in the consecutive plane it
                    % results as an X-crossing. Consequently, in the merged
                    % graph, the two previous branch nodes result connected
                    % to each other. The X-crossing results with degree
                    % 5 and the regular BP with degree 3.

                    % Two options:
                    % Option 1: turn 5-deg node into a 3-deg node by
                    % removing edges between endnodes and 5-deg node.
                    % Option 2: remove edges between endnodes and 3-deg
                    % node and remove edge between 3-deg node and 5-deg
                    % node (i.e, remove all edges of 3-deg node).

                    % Option 1 is the only option that fits with the rest.
                    for iii=1:numel(paths)
                        path = paths{iii};
                        pathDegs = degree(G_merged, path);
                        if any(pathDegs>4)
                            toRemove = path;
                        end
                    end
                end


            end

            for iii=1:numel(toRemove)-1
                G_merged = rmedge(G_merged, toRemove(iii), toRemove(iii+1));
            end

        end


        % Remove tails
        for ii=1:numel(identifiedTails)
            current_tail = identifiedTails{ii};
            current_remainingDiffDegNode = not_in_a_cycle(ii);
            path_focus_values = zeros(numel(current_tail), 1);
            for iii=1:numel(current_tail)
                current_sub_tail = current_tail{iii};
                current_sub_tail = setdiff(current_sub_tail,...
                    current_remainingDiffDegNode);
                idx = findnode(G_test, current_sub_tail);
                path_focus_values(iii) = mean(G_test.Nodes.Node_FM_value(idx));
            end
            if ~isempty(idx) && ~isempty(current_tail)
                [~, idx] = min(path_focus_values);
                toRemove = setdiff(current_tail{idx},current_remainingDiffDegNode);
                G_merged = rmnode(G_merged, toRemove);
            end
        end


        % Remove nodes with degree 0.
        deg0 = G_merged.Nodes.Name(degree(G_merged)==0);
        G_merged = rmnode(G_merged, deg0);


        if debug_mode
            x_coords = G_merged.Nodes.Coordinates(:, 2);
            y_coords = G_merged.Nodes.Coordinates(:, 1);

            % The changed nodes are nodes one should pay attention to, they
            % do not necessarily indicate an errror.

            figure,
            imshow(Img)
            hold on
            p=plot(G_merged, 'XData', x_coords, 'YData', y_coords, 'LineStyle',...
                '-','NodeLabel',{},...
                'EdgeColor', 'w', 'NodeColor', 'r', 'LineWidth', 2, 'MarkerSize', 3);
            

 
            degreeDifference = degree(G_merged) - G_merged.Nodes.OriginalDegree;

            changedNodes = [];

            if sum(abs(degreeDifference)) == 0
                disp("Graph simplification finished")
            else
                changedNodes = G_merged.Nodes.Name(degreeDifference~=0);
                disp("Changed Nodes:")
                disp(changedNodes)
            end

            if ~isempty(changedNodes)
                highlight(p, changedNodes, 'NodeColor', 'r')
            end
            title('Current best graph')
        end

        G_best_so_far = G_merged;

    end

end

end
