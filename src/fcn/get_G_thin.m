function [G_thin, G_minimal] = get_G_thin(sk_structure, densityFactor)

% Returns the representative graph G_thin and the minimal homeomorphic
% graph G_minimal.
% Input:   
%       sk_structure, [N x M] binary skeleton matrix.
%       densityFactor, float 

% Get node coordinates
[I, J] = find(sk_structure);

% Create graph object
G = graph;

% Add nodes to the graph
G = addnode(G, numel(I));

for i = 1:numel(I)
    for j = i+1:numel(I)
        if max(abs([I(i)-I(j), J(i)-J(j)])) <= 1
            G = addedge(G, i, j, hypot(I(i)-I(j), J(i)-J(j)));
        end
    end
end

SG = G;

[cycles, edgecycles] = allcycles(SG, "MaxCycleLength", 3);

edges2remove = [];

for i=1:numel(cycles)
    cycle = cycles{i};
    edgecycle = edgecycles{i};
    if (numel(cycle) == 3) || sum(SG.Edges.Weight(edgecycle) == (2+sqrt(2)))
        edge2remove = edgecycle(SG.Edges.Weight(edgecycle)==sqrt(2));
        edges2remove = [edges2remove; edge2remove];
    end
end

SG = rmedge(SG, edges2remove);


%% Modify graph to desired structure
G_minimal = pruneShortBranches(SG);

G_original = G_minimal;


G_minimal = removeBridgeNodes(G_minimal);

%% Convert indices from char to int

index_tab = table2array(G_minimal.Nodes);
indices = cellfun(@str2num, index_tab);
G_minimal.Nodes.Coordinates = [I(indices), J(indices)];

%% Compute the actual G_minimal
G_minimal = removeBridgeNodes(G_minimal);


%% Create thin graph

subGraphPaths = extractSubgraphPaths(G_original, G_minimal);

G_thin = createThinGraph(G_original, subGraphPaths, densityFactor);

nodes = str2double(G_thin.Nodes.Name);

%% Define points of interests

X_pts = I(nodes);

Y_pts = J(nodes);

test_points = [X_pts, Y_pts];

G_thin.Nodes.Coordinates = test_points;

G_thin.Nodes.Z_planes = ones(size(test_points, 1), 1);

G_thin.Nodes.NodeColors = ones(size(test_points, 1), 1);



end