clear all
rng('default')

addpath(genpath("fcn\"));

warning('off') % Suppress warnings


data_path_prefix = "..\data\Image_Z_stacks_fungi\";

sample_list = ["sample15\", "sample23\"];


sample = 1; % Sample 15 is the sample used to create figure 7

folder = sample_list(sample);
files = dir(fullfile(strcat(data_path_prefix, folder), '*.tif'));

% Parameters
filter_sigma = 5.5; % Depends on length scales
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

% Plane numbers to consider
selectedPlane = [5, 10, 15, 20, 25, 30];

planeNumber = selectedPlane;


% Global variables
G_thin_list = {};
allSporeMasks = zeros(rows, columns, numel(selectedPlane));

%% Run detection for all selected planes

for plane=1:numel(selectedPlane)
    
    % Load image
    img_path = strcat(data_path_prefix, folder,...
        files(selectedPlane(plane)).name);
    current_image_color = imread(img_path);

    % Size computations
    [current_image, mm_per_pixel, position] =...
        scale_bar_removal_color(current_image_color);

    % Run filament detector in the current focal plane
    [sk_structure, ridgeFilt] = mycelium_detection(current_image,...
        polarity, filter_sigma, theta, minimum_branch_length);

    % Remove artifacts introduced by spores
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


    % Remove image borders
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
    
    % Remove small components
    sk_structure = bwareaopen(sk_structure, 10);
    
    % Transform into representative/thin graph
    [G_thin, G_minimal] = get_G_thin(sk_structure, densityFactor);
    
    % Store the graph
    G_thin_list{end+1} = G_thin;


end


%% Plot all the graphs
color_list = ["[0 0.4470 0.7410]", "[0.8500 0.3250 0.0980]",...
    "[0.9290 0.6940 0.1250]", "[0.4940 0.1840 0.5560]", ...
    "[0.4660 0.6740 0.1880]", "[0.3010 0.7450 0.9330]", ...
    "[0.6350 0.0780 0.1840]", "[0 0.6 0]", "[0 0 1]", ...
    "[1 0 0]", "[0 1 1]", "[1 1 1]"];

lgd_list = [];
lgd_str = [];

img_path = strcat(data_path_prefix, folder, files(15).name);
Img = imread(img_path);


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

x_coords = G_updated.Nodes.Coordinates(:, 2);
y_coords = G_updated.Nodes.Coordinates(:, 1);
figure,
imshow(Img)
hold on
plot(G_updated,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b')
title("Correction for spore occlusions")

%% Identify overlaps

G_input = G_updated;
G_minimal = removeBridgeNodes(G_updated);

[nodeIDCrossingPairs, crossings] = identify_overlaps(...
    current_image, G_input, G_minimal, filter_sigma,...
    side_margin, overlap_displacement_tolerance);

% Plot

color_list = ["[0 0.4470 0.7410]", "[0.8500 0.3250 0.0980]",...
    "[0.9290 0.6940 0.1250]", "[0.4940 0.1840 0.5560]", ...
    "[0.4660 0.6740 0.1880]", "[0.3010 0.7450 0.9330]", ...
    "[0.6350 0.0780 0.1840]", "[0 0.6 0]", "[0 0 1]", ...
    "[1 0 0]", "[0 1 1]", "[1 1 1]"];

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
    plot(crossings_tmp(:, 2), crossings_tmp(:, 1), 'LineWidth',...
        3, 'Color', color_list(mod(i, numel(color_list))+1))
end
title("Identified overlaps")


%% Depth interpolation
considered_span = 3;

G_input = depthInterpolation(G_input,...
    reliability_threshold, considered_span, imageMetaData);
%% Match overlaps

G_updated = match_overlaps(G_input, G_minimal, nodeIDCrossingPairs);

%% Visualize
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
    'EdgeColor', 'k', 'NodeCData', colors, 'MarkerSize', 4,...
    'LineWidth', 2, 'EdgeColor', 'w')
colormap('jet');
colorbar;
title("2.5D Representation")

%% 3D Plot

x_coords = G_updated.Nodes.Coordinates(:, 2)*mm_per_pixel;
y_coords = (size(Img, 2)-G_updated.Nodes.Coordinates(:, 1))*mm_per_pixel;
z_coords =  0.02*(max(planeNumber)-G_updated.Nodes.Z_coordinates);

figure,
plot(G_updated, 'XData', x_coords, 'YData', y_coords,...
    'ZData', z_coords, 'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b',...
    'MarkerSize', 4, 'LineWidth', 2)


axis equal
title("3D Representation")


xlabel("[mm]", 'FontSize', 12)
ylabel("[mm]", 'FontSize', 12)
zlabel("[mm]", 'FontSize', 12)

depthRange = max(z_coords)-min(z_coords);

disp(strcat("Depth range: ", string(depthRange), "mm"))






