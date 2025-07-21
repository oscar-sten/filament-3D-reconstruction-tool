% Printed filament reconstruction
% Script for reproducing results in Figure (6) & Table (3)


clear all
rng('default')

addpath(genpath("fcn\"));

warning('off') % Suppress warnings


prototype_number = 2; % Select 1 or 2.

if prototype_number == 1
    data_path_prefix = "..\data\Prototype1\";
    filter_sigma = 18;
    selectedPlane = ["04", "05", "06", "07", "08", "09", "10",...
        "11", "12", "13", "14", "15"];
    planeNumber = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

elseif prototype_number == 2
    data_path_prefix = "..\data\Prototype2\";
    filter_sigma = 30;
    selectedPlane = ["01", "02", "03", "04", "05", "06", "07",...
        "08", "09", "10", "11", "12"];
    planeNumber = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
else
    error("Invalid prototype number")
end



% Common parameters
costUnmatched = 20;
theta = [0:15:360];
polarity = 'dark';
patch_range = [30 30];
minimum_branch_length = 50;
densityFactor = (1/3);
mm_per_pixel = 19/1024;
file_suffix = ".tif";
FM_type = "LAPV";
reliability_threshold = 20; % Reliability threshold on dB for Pertuz's
% criterion.

pruneAngleThreshold = 45;
overlap_displacement_tolerance = 0; % Displacement tolerance in number
% of pixels in each direction.
side_margin = 10;


% Example image
file_prefix = "single_Z-Series_Z";
folder = "";
files = dir(fullfile(strcat(data_path_prefix, folder), '*.tif'));
Img = imread(strcat(data_path_prefix, "single_Z-Series_Z05.tif"));


G_thin_list = {};

%% Run detection



for plane=1:numel(selectedPlane)
    
    % Load image from current focal plane
    img_path = strcat(data_path_prefix, file_prefix,...
        string(selectedPlane(plane)), file_suffix);
    current_image_color = imread(img_path);
    current_image = current_image_color;
    
    % Filament detection
    [sk_structure, ridgeFilt] = mycelium_detection(current_image,...
        polarity, filter_sigma, theta, minimum_branch_length);

    % Transform into representative/thin graph
    [G_thin, G_minimal] = get_G_thin(sk_structure, densityFactor);
    
    % Store the graph
    G_thin_list{end+1} = G_thin;

end


%% Plot
color_list = ["[0 0.4470 0.7410]", "[0.8500 0.3250 0.0980]",...
    "[0.9290 0.6940 0.1250]", "[0.4940 0.1840 0.5560]", ...
    "[0.4660 0.6740 0.1880]", "[0.3010 0.7450 0.9330]", ...
    "[0.6350 0.0780 0.1840]", "[0 0.6 0]", "[0 0 1]",...
    "[0 0 0]", "[1 0 0]", "[0 1 1]", "[1 1 1]"];

lgd_list = [];
lgd_str = [];

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


%% Identify overlaps

G_input = G_merged;
G_minimal = removeBridgeNodes(G_merged);

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

if size(Img, 3)<3
    Img = cat(3, Img, Img, Img); % Convert to rgb image for compatibility with colormap
end

figure,
imshow(Img)
hold on
plot(G_updated, 'XData', x_coords, 'YData', y_coords, 'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'k', 'NodeCData', colors)
colormap('jet');
colorbar;

title("2.5D Representation") 

%% 3D Plot

x_coords = G_updated.Nodes.Coordinates(:, 2)*mm_per_pixel;
y_coords = 19-G_updated.Nodes.Coordinates(:, 1)*mm_per_pixel;
z_coords =  0.5*(max(planeNumber)-G_updated.Nodes.Z_coordinates);

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

depthRange = max(z_coords)-min(z_coords);

disp(strcat("Depth range: ", string(depthRange), "mm"))

%% Ground truth model


if prototype_number == 1
    % Image coords line 1: bottom: [633, 851], top: [727, 227];
    % Here we calculate the leaning of the prototype in the image with respect
    % to the yz-plane. The angle is about 8.57 degrees.
    
    % A "base-depth", keep in mind that we can only estimate relative
    % depth, however, this is sufficient for estimating 3D shape.
    reference_depth = 2.75;

    bottom = [633, 851];
    top = [727, 227];

    % Parallel beam right

    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line1 = zeros(n, 2);

    % Calculate the points
    for i = 1:n
        ref_line1(i, :) = point1 + t(i) * (point2 - point1);
    end

    depth = reference_depth*ones(size(ref_line1, 1), 1);

    ref_line1 = [ref_line1, depth];

    idx_line2 = [23, 24, 25, 10, 26, 27, 28, 29, 31, 32, 33];
    coords_line2 = [x_coords(idx_line2), y_coords(idx_line2),...
        z_coords(idx_line2)];


    [~, D] = dsearchn(ref_line1, coords_line2);


    parallel_line1_dist_RMS = rms(D);

    parallel_line1_dist_sigma = std(D);


    % Other parallel beam

    bottom = [388, 181];
    top = [291, 792];

    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line2 = zeros(n, 2);

    % Calculate the points
    for i = 1:n
        ref_line2(i, :) = point1 + t(i) * (point2 - point1);
    end

    depth = reference_depth*ones(size(ref_line2, 1), 1);

    ref_line2 = [ref_line2, depth];

    idx_line1 = [3, 2, 1, 4, 5, 11, 12, 13, 14, 15, 16];
    coords_line1 = [x_coords(idx_line1), y_coords(idx_line1),...
        z_coords(idx_line1)];


    [~, D] = dsearchn(ref_line2, coords_line1);


    parallel_line2_dist_RMS = rms(D);

    parallel_line2_dist_sigma = std(D);

    % Bridge beam
    idx_bridge = [1, 6, 7, 8, 9, 10];
    coords_bridge = [x_coords(idx_bridge), y_coords(idx_bridge),...
        z_coords(idx_bridge)];

    bottom = [306, 686];
    top = [651, 740];

    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line3 = zeros(n, 2);

    % Calculate the points
    for i = 1:n
        ref_line3(i, :) = point1 + t(i) * (point2 - point1);
    end

    depth = reference_depth*ones(size(ref_line3, 1), 1);

    ref_line3 = [ref_line3, depth];


    [~, D] = dsearchn(ref_line3, coords_bridge);


    bridge_line_dist_RMS = rms(D);

    bridge_line_dist_sigma = std(D);


    % Diagonal beam
    idx = [5, 17, 19, 20, 21, 22, 34, 35, 36];

    coords_diag = [x_coords(idx), y_coords(idx), z_coords(idx)];


    bottom = [321, 599];
    top = [827, 292];

    diff = abs(top-bottom);

    leaning_angle = rad2deg(atan(diff(1)/diff(2)));


    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line4 = zeros(n, 2);

    % Calculate the points
    for i = 1:n
        ref_line4(i, :) = point1 + t(i) * (point2 - point1);
    end

    depth_start = reference_depth;
    angle_xy = 15;

    angle_xy_rad = deg2rad(angle_xy);

    depth = depth_start+...
        (hypot(ref_line4(:, 1)-ref_line4(1, 1),...
        ref_line4(:, 2)-ref_line4(1, 2))*tan(angle_xy_rad));

    ref_line4 = [ref_line4, depth];

    [~, D] = dsearchn(ref_line4, coords_diag);


    diag_beam_dist_RMS = rms(D);

    diag_beam_dist_sigma = std(D);


    all_lines = [ref_line1; ref_line2; ref_line3; ref_line4];
    coords_all = [x_coords, y_coords, z_coords];

    [~,D] = dsearchn(all_lines, coords_all);


    all_dist_RMS = rms(D);

    all_dist_sigma = std(D);

    % Visualize ground truth model
    figure,
    plot3(ref_line1(:, 1), ref_line1(:, 2), ref_line1(:, 3),...
        'b', 'LineWidth', 2)
    hold on
    plot3(ref_line2(:, 1), ref_line2(:, 2), ref_line2(:, 3),...
        'b', 'LineWidth', 2)

    plot3(ref_line3(:, 1), ref_line3(:, 2), ref_line3(:, 3),...
        'b', 'LineWidth', 2)

    plot3(ref_line4(:, 1), ref_line4(:, 2), ref_line4(:, 3),...
        'b', 'LineWidth', 2)
    axis equal

    title("Ground truth model")


    % Visualize reconsturction and ground truth model together
    figure,
    plot3(ref_line1(:, 1), ref_line1(:, 2), ref_line1(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)
    hold on
    plot3(ref_line2(:, 1), ref_line2(:, 2), ref_line2(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)

    plot3(ref_line3(:, 1), ref_line3(:, 2), ref_line3(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)

    plot3(ref_line4(:, 1), ref_line4(:, 2), ref_line4(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)

    x_coords = G_updated.Nodes.Coordinates(:, 2)*mm_per_pixel;
    y_coords = 19-G_updated.Nodes.Coordinates(:, 1)*mm_per_pixel;
    z_coords =  0.5*(max(planeNumber)-G_updated.Nodes.Z_coordinates);

    plot(G_updated, 'XData', x_coords, 'YData', y_coords,...
        'ZData', z_coords, 'LineStyle',...
        '-','NodeLabel',{},...
        'EdgeColor', 'b', 'NodeColor', 'r', 'LineWidth', 1.5, 'MarkerSize', 6)
    axis equal

    legend('Ground truth','' ,'', '', 'Node positions', 'FontSize', 16)

    title("Ground truth model and 3D reconstruction together")


    % Print results
    disp(strcat("------------------ Reconstuction error for prototype #",...
        string(prototype_number), " ----------------"))
    disp(strcat("Filament 1 (Parallell beam on left hand side): RMSE = ", ...
        string(round(parallel_line2_dist_RMS, 2)), " (std = ",...
        string(round(parallel_line2_dist_sigma, 2)), ")"))
    disp(strcat("Filament 2 (Bridge bream): RMSE = ", ...
        string(round(bridge_line_dist_RMS, 2)), " (std = ",...
        string(round(bridge_line_dist_sigma, 2)), ")"))
    disp(strcat("Filament 3 (Parallell beam on right hand side): RMSE = ", ...
        string(round(parallel_line1_dist_RMS, 2)), " (std = ",...
        string(round(parallel_line1_dist_sigma, 2)), ")"))
    disp(strcat("Filament 4 (Diagonal beam): RMSE = ", ...
        string(round(diag_beam_dist_RMS, 2)), " (std = ",...
        string(round(diag_beam_dist_sigma, 2)), ")"))
    disp(strcat("Total: RMSE = ", ...
        string(round(all_dist_RMS, 2)), " (std = ",...
        string(round(all_dist_sigma, 2)), ")"))

elseif prototype_number == 2
    
    % A "base-depth", keep in mind that we can only estimate relative
    % depth, however, this is sufficient for estimating 3D shape.
    reference_depth = 2.25;

    % Vertical beam

    bottom = [475, 1012];
    top = [600, 9];


    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line1 = zeros(n, 2);

    % Calculate the points
    for i = 1:n
        ref_line1(i, :) = point1 + t(i) * (point2 - point1);
    end

    bending = 1; % The filament is slightly bent and we model it as an arch

    if ~bending
        depth = reference_depth*ones(size(ref_line1, 1), 1);
    else
        depth_start = reference_depth;
        top_intersection = [582, 171];
        bottom_intersection = [487, 936];
        bottom_intersection = [bottom_intersection(1)*mm_per_pixel,...
            19-bottom_intersection(2)*mm_per_pixel];
        top_intersection = [top_intersection(1)*mm_per_pixel,...
            19-top_intersection(2)*mm_per_pixel];
        base = norm(top_intersection-bottom_intersection);
        height = 0.4;
        radius = (((base/2)^2)+(height^2))/(2*height);
        base_line = zeros(n, 1);

        % Calculate the points
        for i = 1:n
            base_line(i) = norm(point1) + t(i) * norm(point2 - point1);
        end

        circle_center_x = norm(bottom_intersection)+(base/2)+1;

        circle_center_y = depth_start+height-radius;

        base_line_zerocentered = base_line-circle_center_x;


        depth = -(sqrt((radius^2) - (base_line_zerocentered.^2))...
            -radius+height-depth_start);

    end

    ref_line1 = [ref_line1, depth];

    idx_line1 = findnode(G_updated,...
        string([1, 2, 3, 7, 8, 9, 10, 21, 22, 6, 32, 33, 34, 31, 35, 36, 37]'));
    coords_line1 = [x_coords(idx_line1), y_coords(idx_line1),...
        z_coords(idx_line1)];


    [~, D] = dsearchn(ref_line1, coords_line1);


    vertical_beam_dist_RMS = rms(D);

    vertical_beam_dist_sigma = std(D);


    % Parallel beam left

    bottom = [153, 615];
    top = [582, 171];


    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line2 = zeros(n, 2);

    % Calculate the points
    for i = 1:n
        ref_line2(i, :) = point1 + t(i) * (point2 - point1);
    end

    depth = reference_depth*ones(size(ref_line2, 1), 1);

    ref_line2 = [ref_line2, depth];

    idx_line2 = findnode(G_updated,...
        string([26, 27, 28, 29, 30, 31]'));
    coords_line2 = [x_coords(idx_line2), y_coords(idx_line2),...
        z_coords(idx_line2)];


    [~,D] = dsearchn(ref_line2, coords_line2);


    parallel_beam_left_dist_RMS = rms(D);

    parallel_beam_left_dist_sigma = std(D);


    % Parallel beam right

    bottom = [487, 936];
    top = [912, 467];


    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line3 = zeros(n, 2);

    % Calculate the points
    for i = 1:n
        ref_line3(i, :) = point1 + t(i) * (point2 - point1);
    end

    depth = reference_depth*ones(size(ref_line3, 1), 1);

    ref_line3 = [ref_line3, depth];

    idx_line3 = findnode(G_updated,...
        string([3, 12, 13, 14, 15, 16]'));
    coords_line3 = [x_coords(idx_line3), y_coords(idx_line3),...
        z_coords(idx_line3)];


    [~,D] = dsearchn(ref_line3, coords_line3);


    parallel_beam_right_dist_RMS = rms(D);

    parallel_beam_right_dist_sigma = std(D);

    % Arch

    bottom = [153, 615];
    top = [912, 467];


    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line4 = zeros(n, 2);

    base_line = zeros(n, 1);

    % Calculate the points
    for i = 1:n
        ref_line4(i, :) = point1 + t(i) * (point2 - point1);
        base_line(i) = norm(point1) + t(i) * norm(point2 - point1);
    end

    depth_start = reference_depth;

    % Arch params
    base = 14.42;
    height = 2.65;


    radius = (((base/2)^2)+(height^2))/(2*height);

    circle_center_x = norm(point1)+(base/2);

    circle_center_y = depth_start+height-radius;

    base_line_zerocentered = base_line-circle_center_x;


    depth = sqrt((radius^2) - (base_line_zerocentered.^2))...
        -radius+height+depth_start;


    ref_line4 = [ref_line4, depth];

    idx_line4 = findnode(G_updated,...
        string([16, 17, 18, 19, 11, 20, 23, 24, 25, 26]'));
    coords_line4 = [x_coords(idx_line4), y_coords(idx_line4),...
        z_coords(idx_line4)];


    [~, D] = dsearchn(ref_line4, coords_line4);


    arch_dist_RMS = rms(D);

    arch_dist_sigma = std(D);

    % Visualize ground truth model
    figure,
    plot3(ref_line1(:, 1), ref_line1(:, 2), ref_line1(:, 3),...
        'b', 'LineWidth', 2)
    hold on
    plot3(ref_line2(:, 1), ref_line2(:, 2), ref_line2(:, 3),...
        'b', 'LineWidth', 2)

    plot3(ref_line3(:, 1), ref_line3(:, 2), ref_line3(:, 3),...
        'b', 'LineWidth', 2)

    plot3(ref_line4(:, 1), ref_line4(:, 2), ref_line4(:, 3),...
        'b', 'LineWidth', 2)
    axis equal

    title("Ground truth model")


    all_lines = [ref_line1; ref_line2; ref_line3; ref_line4];
    excluded_nodes = findnode(G_updated, string([4, 5]') );
    included_nodes = setdiff(1:size(x_coords, 1), excluded_nodes);
    coords_all = [x_coords(included_nodes, :),...
        y_coords(included_nodes, :), z_coords(included_nodes, :)];

    [K,D] = dsearchn(all_lines, coords_all);


    all_dist_RMS = rms(D);

    all_dist_sigma = std(D);

    % Visualis GT comparison
    figure,
    plot3(ref_line1(:, 1), ref_line1(:, 2), ref_line1(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)
    hold on
    plot3(ref_line2(:, 1), ref_line2(:, 2), ref_line2(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)

    plot3(ref_line3(:, 1), ref_line3(:, 2), ref_line3(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)

    plot3(ref_line4(:, 1), ref_line4(:, 2), ref_line4(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)


    x_coords = G_updated.Nodes.Coordinates(:, 2)*mm_per_pixel;
    y_coords = 19-G_updated.Nodes.Coordinates(:, 1)*mm_per_pixel;
    z_coords =  0.5*(max(planeNumber)-G_updated.Nodes.Z_coordinates);

    plot(G_updated, 'XData', x_coords, 'YData', y_coords,...
        'ZData', z_coords, 'LineStyle',...
        '-','NodeLabel',{},...
        'EdgeColor', 'b', 'NodeColor', 'r', 'LineWidth', 1.5, 'MarkerSize', 6)
    axis equal

    legend('Ground truth','' ,'', '', '3D Reconsruction', 'FontSize', 16)

    title("Ground truth model and 3D reconstruction together")

    % Print results
    disp(strcat("------------------ Reconstuction error for prototype #",...
        string(prototype_number), " ----------------"))
    disp(strcat("Filament 1 (Parallell beam on left hand side): RMSE = ", ...
        string(round(parallel_beam_left_dist_RMS, 2)), " (std = ",...
        string(round(parallel_beam_left_dist_sigma, 2)), ")"))
    disp(strcat("Filament 2 (Parallell beam on right hand side): RMSE = ", ...
        string(round(parallel_beam_right_dist_RMS, 2)), " (std = ",...
        string(round(parallel_beam_right_dist_sigma, 2)), ")"))
    disp(strcat("Filament 3 (Vertical beam): RMSE = ", ...
        string(round(vertical_beam_dist_RMS, 2)), " (std = ",...
        string(round(vertical_beam_dist_sigma, 2)), ")"))
    disp(strcat("Filament 4 (Arch): RMSE = ", ...
        string(round(arch_dist_RMS, 2)), " (std = ",...
        string(round(arch_dist_sigma, 2)), ")"))
    disp(strcat("Average: RMSE = ", ...
        string(round(all_dist_RMS, 2)), " (std = ",...
        string(round(all_dist_sigma, 2)), ")"))

else
    error("Invalid prototype number")
end




