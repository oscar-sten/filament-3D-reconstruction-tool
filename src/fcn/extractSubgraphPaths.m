function subGraphPaths = extractSubgraphPaths(G_original, G_minimal)

subGraphPaths = {}; % Cell with all line segements

G_original_cpy = G_original;

test = G_minimal.Edges.EndNodes;

numArray = cellfun(@str2num, test);

% Find unique rows in numArray
[uniqueRows, idx, ~] = unique(numArray, 'rows', 'stable');

% Find non-unique rows in numArray
nonUniqueRows = setdiff(1:size(numArray, 1), idx);

nonUniqueEdges = numArray(nonUniqueRows, :);

edges_base = uniqueRows;

[n_node_pairs, ~] = size(edges_base);

for i=1:n_node_pairs
    node_pair = edges_base(i, :);
    if ~ismember(node_pair, nonUniqueEdges, 'rows')
        path = shortestpath(G_original, string(node_pair(1)),...
            string(node_pair(2)));
        degree3 = (degree(G_original, path)==3);
        if sum(degree3)>2
            path = setdiff(path', string(node_pair'));
            G_original_cpy = rmnode(G_original_cpy, path);
            secondPath = shortestpath(G_original_cpy, string(node_pair(1)), ...
                string(node_pair(2)));
            path = secondPath;
        end
        subGraphPaths{end+1} = path;
    else
        [~, ~, idx] = intersect(node_pair, nonUniqueEdges, 'rows');
        n_equal_rows = numel(idx);

        for j=1:1+n_equal_rows
            
%             % Current
%             path = shortestpath(G_original_cpy, string(node_pair(1)), ...
%                 string(node_pair(2)));
% 
%             degree3 = (degree(G_original, path)==3);
%             if sum(degree3)==2
%                 subGraphPaths{end+1} = path;
%                 path = setdiff(path', string(node_pair'));
%                 G_original_cpy = rmnode(G_original_cpy, path);
%             else
%                 path = setdiff(path', string(node_pair'));
%                 G_original_cpy = rmnode(G_original_cpy, path);
%                 path = shortestpath(G_original_cpy, string(node_pair(1)), ...
%                     string(node_pair(2)));
%                 subGraphPaths{end+1} = path;
%             end
            
            % Experiamental       
            try
                path = shortestpath(G_original_cpy, string(node_pair(1)), ...
                    string(node_pair(2)));
            catch
                path = shortestpath(G_original, string(node_pair(1)), ...
                    string(node_pair(2)));
            end
            % Temporary solution:
            degree3 = (degree(G_original, path)>=3);
            if sum(degree3)==2
                subGraphPaths{end+1} = path;
                path = setdiff(path', string(node_pair'));
                G_original_cpy = rmnode(G_original_cpy, path);
            else
                path = setdiff(path', string(node_pair'));
                G_original_cpy = rmnode(G_original_cpy, path);
                path = shortestpath(G_original_cpy, string(node_pair(1)), ...
                    string(node_pair(2)));
                subGraphPaths{end+1} = path;
            end



        end
    end
end



end

