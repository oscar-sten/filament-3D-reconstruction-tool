
% Code for reproducing Figures 8 and 9.
% The change the variable "sample_number".
% sample_number = 1 -> Figure 7
% sample_number = 2 -> Figure 8 & 9
% sample_number = 3 -> Figure S9


clear all
rng('default')

addpath(genpath("fcn\"));

warning('off') % Suppress warnings


data_path_prefix = "..\data\Image_Z_stacks_fungi\";

sample_number = 1; % Set to 1, 2, or 3.

if sample_number == 1
    folder = "sample15\";
    filter_sigma = 5.5; % Depends on length scales
    
    % Plane numbers to consider
    selectedPlane = [5, 10, 15, 20, 25, 30];
    planeNumber = selectedPlane;

elseif sample_number == 2
    folder = "sample23\";
    filter_sigma = 4.5; % Depends on length scales

    % Plane numbers to consider
    selectedPlane = [5, 10, 15, 20, 25];
    planeNumber = selectedPlane;

elseif sample_number == 3
    folder = "sample71\";
    filter_sigma = 4.5; % Depends on length scales

    selectedPlane = [10, 15, 20, 25, 28, 29];
    planeNumber = selectedPlane;

end

% folder = sample_list(sample);
files = dir(fullfile(strcat(data_path_prefix, folder), '*.tif'));

% Parameters (same for all)
theta = [0:15:360]; % Filtering resolution, not changed.
polarity = 'dark';  % Polarity of filaments, usually dark
patch_range = [10 10]; % Sets window size. [x y] indicates how many 
% pixlels in each direction from the node coordinates to considider.
minimum_branch_length = 20; % Depends on length scales
costUnmatched = 3; % Depends on length scales
densityFactor = 1;
pruneAngleThreshold = 45;

FM_type = "TENV"; % 4 letter combination used as argument for the function
% 'fmeasure.m' originaly written by Pertuz et al. (2013),
% reproduced in fcn/external under Open Source BSD licence.
reliability_threshold = 25; % Reliability threshold on dB for Pertuz's
% criterion.
overlap_displacement_tolerance = 0; % Displacement tolerance in number
% of pixels in each direction.
side_margin = 10;

file_suffix = ".tif";

% Scale ...
img_path = strcat(data_path_prefix, folder, files(1).name);
current_image_color = imread(img_path);

% Size computations
[current_image, mm_per_pixel] = scale_bar_removal_color(...
    current_image_color);
spore_diameter_range_mm = [0.025 0.14];

lower_lim_pix = (spore_diameter_range_mm(1)/2)/mm_per_pixel;
upper_lim_pix = (spore_diameter_range_mm(2)/2)/mm_per_pixel;

lower_lim = 6;


% scale_factor = 4;
scale_factor = lower_lim_pix/lower_lim;

% uppper_lim = 30;
uppper_lim = ceil(upper_lim_pix/scale_factor);

patch_size_mm = 0.05;
patch_size = round(patch_size_mm/mm_per_pixel);

% CHT parameters
sensitivity = 0.8;
method = "PhaseCode";

[rows, columns] = size(current_image);
sk_stack = zeros(rows, columns, 2);


% Global variables
G_thin_list = {};
allSporeMasks = zeros(rows, columns, numel(selectedPlane));


%% Run detection



for plane=1:numel(selectedPlane)

    img_path = strcat(data_path_prefix, folder,...
        files(selectedPlane(plane)).name);
    current_image_color = imread(img_path);

    % Size computations
    [current_image, mm_per_pixel, position] =...
        scale_bar_removal_color(current_image_color);


    [sk_structure, ridgeFilt] = mycelium_detection(current_image,...
        polarity, filter_sigma, theta, minimum_branch_length);


    current_image_rescaled = imresize(current_image, 1/scale_factor);
    [centers, radii, metric] = imfindcircles(current_image_rescaled,...
        [lower_lim uppper_lim], "ObjectPolarity", "dark",...
        "Sensitivity", sensitivity, "Method", method);

    centers = scale_factor*centers;
    radii = scale_factor*radii;

    [binSporeMask, idSporeMask] = genSporeMask2(centers, radii, rows,...
        columns, sk_structure);
    
    allSporeMasks(:, :, plane) = idSporeMask;

    sk_structure = sk_structure & ~binSporeMask;

    frame = zeros(size(current_image));
    n = 10;
    frame(:, 1:n) = 1;
    frame(:, columns-n:columns) = 1;
    frame(1:n, :) = 1;
    frame(rows-n:rows, :) = 1;

    scale_bar_effect = zeros(size(current_image));
    scale_bar_effect((position(1)-20):(position(2)+20),...
        (position(3)-20):(position(4)+20)) = 1;

    sk_structure = sk_structure & ~frame;
    sk_structure = sk_structure & ~scale_bar_effect;

    sk_structure = bwareaopen(sk_structure, 10);

    % Transform into representative/thin graph
    [G_thin, G_minimal] = get_G_thin(sk_structure, densityFactor);

    G_thin_list{end+1} = G_thin;


end


%% Plot
color_list = ["[0 0.4470 0.7410]", "[0.8500 0.3250 0.0980]",...
    "[0.9290 0.6940 0.1250]", "[0.4940 0.1840 0.5560]", ...
    "[0.4660 0.6740 0.1880]", "[0.3010 0.7450 0.9330]", ...
    "[0.6350 0.0780 0.1840]", "[0 0.6 0]", "[0 0 1]", ...
    "[1 0 0]", "[0 1 1]", "[1 1 1]"];

lgd_list = [];
lgd_str = [];

img_path = strcat(data_path_prefix, folder, files(15).name);
Img = imread(img_path);


% Closest JET-color: [0, 0, 0.6], furthest: [0.6, 0, 0]

figure,
imshow(Img)
hold on
for i=1:numel(G_thin_list)
    current_graph = G_thin_list{i};
    y_coords = current_graph.Nodes.Coordinates(:, 1);
    x_coords = current_graph.Nodes.Coordinates(:, 2);
    lgd_list = [lgd_list; plot(current_graph,'XData', x_coords,'YData',...
        y_coords,'LineStyle',...
        '-','NodeLabel',{},...
        'EdgeColor', 'w', 'NodeColor',...
        color_list(mod(i, numel(color_list))+1))];
    lgd_str{end+1} = strcat("Plane: ", string(selectedPlane(i)));
end

legend(lgd_list, lgd_str)
title('All graphs from the considered focal planes')

%% Merge graphs

imageMetaData = struct;
imageMetaData.dataPathPrefix = data_path_prefix;
imageMetaData.folder = folder;
imageMetaData.files = files;
imageMetaData.selectedPlane = selectedPlane;
imageMetaData.planeNumber = planeNumber;

depthParams = struct;
depthParams.patchRange = patch_range;
depthParams.FMType = FM_type;
depthParams.costUnmatched = costUnmatched;

G_merged = mergeGraphs(G_thin_list, imageMetaData,...
    depthParams, pruneAngleThreshold);


%% Plot final output
x_coords = G_merged.Nodes.Coordinates(:, 2);
y_coords = G_merged.Nodes.Coordinates(:, 1);
figure,
imshow(Img)
hold on
plot(G_merged,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b')
title('Final merged 2D graph')

%% Integrate part from MetroAgriFor

G_updated = G_merged;

[nodeSets, armVectorSet, connectionVectorSet] =...
    computeSets4SporeOcclusionCorrection(G_updated,...
    idSporeMask, allSporeMasks);

%% Plot

x_coords = G_updated.Nodes.Coordinates(:, 2);
y_coords = G_updated.Nodes.Coordinates(:, 1);
figure,
imshow(Img)
hold on
plot(G_updated,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b')
for i=1:numel(nodeSets)
    currentNodeSet = nodeSets{i};
    if numel(currentNodeSet) > 1
        armVectors = armVectorSet{i};
        for j=1:numel(currentNodeSet)
            baseCoords = G_updated.Nodes.Coordinates(currentNodeSet(j), :);
            armVector = armVectors(j, :);
            purlpe_quiver = quiver(baseCoords(2),...
                baseCoords(1),...
                armVector(2), armVector(1),...
                'Color', '[0.4940 0.1840 0.5560]',...
                'LineWidth', 2, 'MaxHeadSize', 5);
        end
        connectionVectors = connectionVectorSet{i};
        for j=1:numel(currentNodeSet)
            baseCoords = G_updated.Nodes.Coordinates(currentNodeSet(j), :);
            for k=1:numel(currentNodeSet)-1
                connectionVector = connectionVectors(j, :, k);
                orange_quiver = quiver(baseCoords(2),...
                    baseCoords(1),...
                    connectionVector(2), connectionVector(1),...
                    'Color', '[0.8500 0.3250 0.0980]',...
                    'LineWidth', 2, 'MaxHeadSize', 5);
            end
        end
        green_star = plot(G_updated.Nodes.Coordinates(currentNodeSet, 2), ...
            G_updated.Nodes.Coordinates(currentNodeSet, 1), '*',...
            'LineWidth', 2, 'Color', "[0, 0.6, 0]");

    else
        coords = G_updated.Nodes.Coordinates(currentNodeSet, :);
        red_star = plot(coords(2), coords(1), 'r*', 'LineWidth', 2);
    end
end

legend([red_star, green_star, purlpe_quiver, orange_quiver], ...
    {"Connection to Spore", "Open Connection", "Arm Vector",...
    "Candidate Connection Vector"}, ...
    'FontSize', 14, 'Location', 'northwest')
title('Integration of MetroAgriFor algorithm')

%% Matching step in MetroAgriFor algorithm

G_updated = sporeOcclusionCorrection(G_updated, ...
    nodeSets, armVectorSet, connectionVectorSet);

%% Plot
x_coords = G_updated.Nodes.Coordinates(:, 2);
y_coords = G_updated.Nodes.Coordinates(:, 1);
figure,
imshow(Img)
hold on
plot(G_updated,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b')
title("Correction for spore occlusions")


% %% In 2D no merging is done, all nodes derive from the same plane.
% G_thin = G_thin_list{3};
% G_thin.Nodes.Z_coordinates = ones(numel(G_thin.Nodes.Z_planes), 1);
% 
% G_thin.Nodes.SourceGraphIdx = ones(numel(G_thin.Nodes.Z_planes), 1);
% G_updated = G_thin;
%% Identify overlaps

G_input = G_updated;
G_minimal = removeBridgeNodes(G_updated);

[nodeIDCrossingPairs, crossings] = identify_overlaps(...
    current_image, G_input, G_minimal, filter_sigma,...
    side_margin, overlap_displacement_tolerance);


%% Plot
x_coords = G_input.Nodes.Coordinates(:, 2);
y_coords = G_input.Nodes.Coordinates(:, 1);

figure,
imshow(Img)
hold on
p=plot(G_input, 'XData', x_coords, 'YData', y_coords, 'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b');

n_test = size(crossings, 1);

for i=1:2:n_test-1
    crossings_tmp = crossings(i:i+1, :);
    plot(crossings_tmp(:, 2), crossings_tmp(:, 1), 'LineWidth', 3)
end

title("Identified overlaps")

%% Depth interpolation
considered_span = 3;

G_input = depthInterpolation(G_input,...
    reliability_threshold, considered_span, imageMetaData);

%% Match overlaps
G_updated = match_overlaps(G_input, G_minimal, nodeIDCrossingPairs);

%% Visualize G_updated
x_coords = G_updated.Nodes.Coordinates(:, 2);
y_coords = G_updated.Nodes.Coordinates(:, 1);

figure,
imshow(Img)
hold on 
plot(G_updated,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b')
title("Graph with overlaps")

%% Smoothing Depths
% This is done after correction for overlaps

G_updated = smoothDepthCoords(G_updated);

%% Plot
x_coords = G_updated.Nodes.Coordinates(:, 2);
y_coords = G_updated.Nodes.Coordinates(:, 1);

colors = G_updated.Nodes.Z_coordinates;
% colors = G_updated.Nodes.Z_planes;

if size(Img, 3)<3
    Img = cat(3, Img, Img, Img); % Convert to rgb image for compatibility with colormap
end

figure,
imshow(Img)
hold on
plot(G_updated, 'XData', x_coords, 'YData', y_coords, 'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'k', 'NodeCData', colors, 'MarkerSize', 2,...
    'LineWidth', 2, 'EdgeColor', 'w')
colormap('jet');
colorbar("Ticks", linspace(min(colors), max(colors), 5),...
    "Ticklabels", {"Closer","","","","Further"},...
    "TickLabelInterpreter", "latex", "FontSize", 14)
title("2.5D Representation")


%% 3D Plot

x_coords = G_updated.Nodes.Coordinates(:, 2)*mm_per_pixel;
y_coords = G_updated.Nodes.Coordinates(:, 1)*mm_per_pixel;
z_coords =  0.02*(max(planeNumber)-G_updated.Nodes.Z_coordinates);

figure,
plot(G_updated, 'XData', x_coords, 'YData', y_coords,...
    'ZData', z_coords, 'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b')


axis equal


xlabel("[mm]")
ylabel("[mm]")
zlabel("[mm]")

title("3D Representation")

%% Divide graph


G_updated2 = G_updated;

% identify deg 2

deg2Nodes = find(degree(G_updated2)==2);


bridgeNodes = G_updated2.Nodes.Name(deg2Nodes);

for i=1:numel(bridgeNodes)
    currentNode = bridgeNodes(i);

    if degree(G_updated2, currentNode)==2
        nghbrs = neighbors(G_updated2, currentNode);
        try
        [~, dist] = shortestpath(G_updated2, nghbrs{1},  nghbrs{2});
        G_updated2 = addedge(G_updated2, nghbrs{1},  nghbrs{2}, dist);
        G_updated2 = rmnode(G_updated2, currentNode);
        catch
        end
    end

end

x_coords = G_updated2.Nodes.Coordinates(:, 2);
y_coords = G_updated2.Nodes.Coordinates(:, 1);

figure,
imshow(Img)
hold on 
plot(G_updated2,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b')

title("Morphological graph with overlaps")

%% Identify subgraphs
subGraphPaths = extractSubgraphPaths(G_updated, G_updated2);


%% Qualitative demo results: length, width, and Aguilar-Trigueros's model


G_edge_lengths = G_updated2;

G_edge_widths = G_updated2;

G_r4 = G_updated2;

for i=1:numel(subGraphPaths)
    sg = subGraphPaths{i};
    sgIdx = findnode(G_updated, sg);
    x_coords = G_updated.Nodes.Coordinates(sgIdx, 1);
    y_coords = G_updated.Nodes.Coordinates(sgIdx, 2);
    coords = [x_coords, y_coords];
    n_coords = size(coords, 1);
    widthData = zeros(n_coords, 1);
    edgeLength = G_updated2.Edges.Weight(i)*mm_per_pixel; % [mm]
    for j=1:n_coords
        currentCoords = coords(j, :);
        patch_range = [15, 15];
        [patch, x_range_1, y_range_1] =...
            extractPatch(current_image, currentCoords, patch_range);

        [mean_width_px, color] = width_measure(patch, 'Otsu');

        mean_width = mean_width_px*mm_per_pixel*1000; % Mean width in microns.

        widthData(j) = mean_width;
    end
    filamentRadius = mean(widthData)/2; % [microns]
    G_edge_lengths.Edges.Weight(i) = edgeLength; % [mm]
    G_edge_widths.Edges.Weight(i) = 2*filamentRadius; %[microns]
    G_r4.Edges.Weight(i) = (1000*edgeLength)/(filamentRadius^4); %[microns*microns^-4]
end



%% Plots



x_coords = G_r4.Nodes.Coordinates(:, 2);
y_coords = G_r4.Nodes.Coordinates(:, 1);


figure,
% subplot(1, 3, 1)
imshow(Img)
hold on 
% h=plot(G_edge_lengths,'XData', x_coords,'YData',y_coords,'LineStyle',...
%     '-','NodeLabel',{},...
%     'EdgeColor', 'r', 'NodeColor', 'b', 'EdgeLabel',...
%     round(G_edge_lengths.Edges.Weight, 2));
h=plot(G_edge_lengths,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
     'NodeColor', 'm', 'EdgeLabel',...
    round(G_edge_lengths.Edges.Weight, 2), 'EdgeCData',...
    G_edge_lengths.Edges.Weight,...
    'LineWidth', 2);
colormap('jet');
colorbar("Ticks", linspace(min(G_edge_lengths.Edges.Weight),...
    max(G_edge_lengths.Edges.Weight),...
    5),"Ticklabels",...
        {"Shorter","","","", "Longer"},...
        "TickLabelInterpreter", "latex", "FontSize", 14);
h.EdgeFontSize=12;
title('Lengths [mm]')


%% dlkö

figure,
% subplot(1, 3, 2)
imshow(Img)
hold on 
h=plot(G_edge_widths,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'b', 'NodeColor', 'm', 'EdgeLabel',...
    round(G_edge_widths.Edges.Weight, 2), 'LineWidth', 0.3*G_edge_widths.Edges.Weight);
h.EdgeFontSize=12;
title('Filament diameter [\mum]')


%% dksl

figure,
% subplot(1, 3, 3)
imshow(Img)
hold on 
h=plot(G_r4,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
   'NodeColor', 'm', 'EdgeLabel',...
    round(G_r4.Edges.Weight, 2), 'EdgeCData', G_r4.Edges.Weight,...
    'LineWidth', 2);
colormap('jet');
colorbar("Ticks", linspace(min(G_r4.Edges.Weight), max(G_r4.Edges.Weight),...
    5),"Ticklabels",...
        {"Lower","","", "","Higher"},...
        "TickLabelInterpreter", "latex", "FontSize", 14)
h.EdgeFontSize=12;
title('Predicted transport resistance [\mu\cdot\mum^{-4}]')
