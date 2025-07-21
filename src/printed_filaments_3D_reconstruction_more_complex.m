% Printed filament reconstruction

% The following code if for reproducing the 3D reconstructions of the
% Prototypes #3-5. These Z-stacks were aquired after the microscopy
% platform had been modified, which is why they have a different scale.
% Because the scale was too large for ridge detection (stuctures resulted
% more like blobs) the detection procedure was slightly modified.
% Nevertheless, the graph merging, overlap detection and depth estimation
% (which are the focus of this paper) are the same.

clear all
rng('default')

addpath(genpath("fcn\"));

warning('off') % Suppress warnings


prototype_number = 3; % Prototype number, choose 3, 4, or 5.


% Specific parameters
if prototype_number == 3
    data_path_prefix = "..\data\Prototype3\";
    file_prefix = "obj4001z";
    filter_sigma = 45;
    strel_radius = 15;
    selectedPlane = ["01", "02", "03", "04", "05", "06", "07", "08",...
        "09", "10"];
    planeNumber = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

elseif prototype_number == 4
    data_path_prefix = "..\data\Prototype4\";
    file_prefix = "obj4002z";
    filter_sigma = 45;
    strel_radius = 15;
    selectedPlane = ["01", "02", "03", "04", "05", "06", "07", "08",...
        "09", "10", "11", "12", "13", "14"];
    planeNumber = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];

elseif prototype_number == 5
    data_path_prefix = "..\data\Prototype5\";
    file_prefix = "obj4005z";
    filter_sigma = 45;
    strel_radius = 5;
    selectedPlane = ["03", "04", "05", "06", "07", "08",...
        "09", "10", "11", "12", "13", "14"];
    planeNumber = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];

end



% Common parameters
costUnmatched = 20;
theta = [0:15:360];
polarity = 'dark';
patch_range = [25 25];
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
side_margin = 20;

folder = "";
files = dir(fullfile(strcat(data_path_prefix, folder), '*.tif'));
Img = imread(strcat(data_path_prefix, file_prefix, "05.tif"));


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
    
    test_alt = imbinarize(mat2gray(ridgeFilt));
    test_alt = imclose(test_alt, strel('disk', strel_radius));
    sk_test_alt = bwskel(test_alt, 'MinBranchLength',...
        minimum_branch_length);

    % Transform into representative/thin graph
    [G_thin, G_minimal] = get_G_thin(sk_test_alt, densityFactor);
    
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

mm_per_pixel = 2/137;

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



%% Model comparison


if prototype_number == 3
    mm_per_pixel = 2/137;

    reference_depth = 0.75;
    
    % Parallel beam left
    bottom = [339, 991];
    top = [218, 102];
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
    idx_line1 = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 19, 20, 21,...
        26, 27]; % 28 outside
    coords_line1 = [x_coords(idx_line1), y_coords(idx_line1),...
        z_coords(idx_line1)];

    [~, D] = dsearchn(ref_line1, coords_line1);

    parallel_line1_dist_RMS = rms(D);

    parallel_line1_dist_sigma = std(D);

    % Other parallel beam

    bottom = [794, 939];
    top = [676, 38];

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

    idx_line2 = [49, 48, 50, 51, 52, 43, 53, 54, 55, 47, 56, 57, 58, 34,...
        59, 60, 61]; 
    coords_line2 = [x_coords(idx_line2), y_coords(idx_line2),...
        z_coords(idx_line2)];

    [~, D] = dsearchn(ref_line2, coords_line2);

    parallel_line2_dist_RMS = rms(D);

    parallel_line2_dist_sigma = std(D);

    % Perpendicular bridge
    bottom = [316, 851];
    top = [773, 791];

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

    idx_line3 = [21, 29, 30, 31, 32, 33, 34]; 
    coords_line3 = [x_coords(idx_line3), y_coords(idx_line3),...
        z_coords(idx_line3)];

    [~, D] = dsearchn(ref_line3, coords_line3);

    horizontal_line_dist_RMS = rms(D);

    horizontal_line_dist_sigma = std(D);


    % Diagonal bridge beam
    bottom = [271, 511];
    top = [688, 109];

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

    depth = reference_depth*ones(size(ref_line4, 1), 1);

    ref_line4 = [ref_line4, depth];

    idx_line4 = [8, 16, 17, 35, 36, 37, 38]; 
    coords_line4 = [x_coords(idx_line4), y_coords(idx_line4),...
        z_coords(idx_line4)];

    [~, D] = dsearchn(ref_line4, coords_line4);

    diagonal_line_dist_RMS = rms(D);

    diagonal_line_dist_sigma = std(D);

    % Small Arch

    bottom = [302, 737];
    top = [716, 333];


    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line5 = zeros(n, 2);

    base_line = zeros(n, 1);

    % Calculate the points
    for i = 1:n
        ref_line5(i, :) = point1 + t(i) * (point2 - point1);
        base_line(i) = norm(point1) + t(i) * norm(point2 - point1);
    end

    depth_start = reference_depth;

    % Arch params
    base = 8.4446;
    height = 1.37;


    radius = (((base/2)^2)+(height^2))/(2*height);

    circle_center_x = norm(point1)+(base/2);

    circle_center_y = depth_start+height-radius;

    base_line_zerocentered = base_line-circle_center_x;


    depth = sqrt((radius^2) - (base_line_zerocentered.^2))...
        -radius+height+depth_start;

    ref_line5 = [ref_line5, depth];

    idx_line5 = [15, 22, 23, 24, 41, 42, 43]; 
    coords_line5 = [x_coords(idx_line5), y_coords(idx_line5),...
        z_coords(idx_line5)];

    [~, D] = dsearchn(ref_line5, coords_line5);

    small_arch_dist_RMS = rms(D);

    small_arch_dist_sigma = std(D);

    % Big Arch

    bottom = [762, 667];
    top = [232, 189];
    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];
    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line6 = zeros(n, 2);

    base_line = zeros(n, 1);

    % Calculate the points
    for i = 1:n
        ref_line6(i, :) = point1 + t(i) * (point2 - point1);
        base_line(i) = norm(point1) + t(i) * norm(point2 - point1);
    end

    depth_start = reference_depth;

    % Arch params
    base = 10.4191;
    height = 3;


    radius = (((base/2)^2)+(height^2))/(2*height);

    circle_center_x = norm(point1)+(base/2);

    circle_center_y = depth_start+height-radius;

    base_line_zerocentered = base_line-circle_center_x;


    depth = sqrt((radius^2) - (base_line_zerocentered.^2))...
        -radius+height+depth_start;

    ref_line6 = [ref_line6, depth];

    idx_line6 = [4, 9, 10, 11, 18, 39, 40, 25, 44, 45, 46, 47]; 
    coords_line6 = [x_coords(idx_line6), y_coords(idx_line6),...
        z_coords(idx_line6)];

    [~, D] = dsearchn(ref_line6, coords_line6);

    big_arch_dist_RMS = rms(D);

    big_arch_dist_sigma = std(D);

    % All 
    all_lines = [ref_line1; ref_line2; ref_line3; ref_line4;...
        ref_line5; ref_line6];
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
    plot3(ref_line5(:, 1), ref_line5(:, 2), ref_line5(:, 3),...
        'b', 'LineWidth', 2)
    plot3(ref_line6(:, 1), ref_line6(:, 2), ref_line6(:, 3),...
        'b', 'LineWidth', 2)
    axis equal

    title("Ground truth model")


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
    plot3(ref_line5(:, 1), ref_line5(:, 2), ref_line5(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)
    plot3(ref_line6(:, 1), ref_line6(:, 2), ref_line6(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)


    x_coords = G_updated.Nodes.Coordinates(:, 2)*mm_per_pixel;
    y_coords = 19-G_updated.Nodes.Coordinates(:, 1)*mm_per_pixel;
    z_coords =  0.5*(max(planeNumber)-G_updated.Nodes.Z_coordinates);

    plot(G_updated, 'XData', x_coords, 'YData', y_coords,...
        'ZData', z_coords, 'LineStyle',...
        '-','NodeLabel',{},...
        'EdgeColor', 'b', 'NodeColor', 'r',...
        'LineWidth', 1.5, 'MarkerSize', 6)
    axis equal

    legend('Ground truth','' ,'', '', '', '', '3D Reconsruction',...
        'FontSize', 16)

    title("Ground truth model and 3D reconstruction together")

    disp(strcat("------------------ Reconstuction error for prototype #",...
        string(prototype_number), " ----------------"))
    disp(strcat("Filament 1 (Parallell beam on left hand side): RMSE = ", ...
        string(round(parallel_line1_dist_RMS, 2)), " (std = ",...
        string(round(parallel_line1_dist_sigma, 2)), ")"))
    disp(strcat("Filament 2 (Parallell beam on right hand side): RMSE = ", ...
        string(round(parallel_line2_dist_RMS, 2)), " (std = ",...
        string(round(parallel_line2_dist_sigma, 2)), ")"))
    disp(strcat("Filament 3 (Bridge bream): RMSE = ", ...
        string(round(horizontal_line_dist_RMS, 2)), " (std = ",...
        string(round(horizontal_line_dist_sigma, 2)), ")"))
    disp(strcat("Filament 4 (Diagonal bream): RMSE = ", ...
        string(round(diagonal_line_dist_RMS, 2)), " (std = ",...
        string(round(diagonal_line_dist_sigma, 2)), ")"))
    disp(strcat("Filament 5 (Small arch): RMSE = ", ...
        string(round(small_arch_dist_RMS, 2)), " (std = ",...
        string(round(small_arch_dist_sigma, 2)), ")"))
    disp(strcat("Filament 6 (Big arch): RMSE = ", ...
        string(round(big_arch_dist_RMS, 2)), " (std = ",...
        string(round(big_arch_dist_sigma, 2)), ")"))
    disp(strcat("Average: RMSE = ", ...
        string(round(all_dist_RMS, 2)), " (std = ",...
        string(round(all_dist_sigma, 2)), ")"))




elseif prototype_number == 4
    mm_per_pixel = 2/137;

    reference_depth = 0.5; % 1.75-1.5;
    
    % Parallel beam left
    bottom = [310, 988];
    top = [281, 53];
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
    idx_line1 = [11, 12, 13, 14, 18, 17, 16, 15, 19, 20, 21, 22, 7, 8, 9]; % 10 outside
    coords_line1 = [x_coords(idx_line1), y_coords(idx_line1),...
        z_coords(idx_line1)];

    [~, D] = dsearchn(ref_line1, coords_line1);

    parallel_line1_dist_RMS = rms(D);

    parallel_line1_dist_sigma = std(D);

    % Other parallel beam

    bottom = [823, 915];
    top = [806, 30];

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

    idx_line2 = [50, 49, 48, 36, 45, 46, 47, 40, 51, 52, 53, 54, 44, 55, 56]; % 57, 58, 59 
    coords_line2 = [x_coords(idx_line2), y_coords(idx_line2),...
        z_coords(idx_line2)];

    [~, D] = dsearchn(ref_line2, coords_line2);

    parallel_line2_dist_RMS = rms(D);

    parallel_line2_dist_sigma = std(D);

        % Perpendicular bridge
    bottom = [289, 212];
    top = [808, 206];

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

    idx_line3 = [14, 27, 28, 29, 30, 31, 32]; 
    coords_line3 = [x_coords(idx_line3), y_coords(idx_line3),...
        z_coords(idx_line3)];

    [~, D] = dsearchn(ref_line3, coords_line3);

    horizontal_line_dist_RMS = rms(D);

    horizontal_line_dist_sigma = std(D);


    % Diagonal beam

    bottom = [811, 328];
    top = [202, 817];

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
    angle_xy = 30;

    angle_xy_rad = deg2rad(angle_xy);

        depth = depth_start+...
            (hypot(ref_line4(:, 1)-ref_line4(1, 1),...
            ref_line4(:, 2)-ref_line4(1, 2))*tan(angle_xy_rad));

    ref_line4 = [ref_line4, depth];

    idx = [40, 39, 38, 37, 35, 34, 33, 6, 5, 4, 3, 2, 1];

    coords_diag = [x_coords(idx), y_coords(idx), z_coords(idx)];


    [~, D] = dsearchn(ref_line4, coords_diag);


    diag_beam_dist_RMS = rms(D);

    diag_beam_dist_sigma = std(D);
    
    % Small Arch

    bottom = [812, 757];
    top = [285, 340];


    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line5 = zeros(n, 2);

    base_line = zeros(n, 1);

    % Calculate the points
    for i = 1:n
        ref_line5(i, :) = point1 + t(i) * (point2 - point1);
        base_line(i) = norm(point1) + t(i) * norm(point2 - point1);
    end
    

    % Arch params
    base = 9.8106;
    height = 1.47;
    
    arch_elevation = 2.5;
    
    depth_start = reference_depth+arch_elevation;

    radius = (((base/2)^2)+(height^2))/(2*height);

    circle_center_x = norm(point1)+(base/2);

    circle_center_y = depth_start+height-radius;

    base_line_zerocentered = base_line-circle_center_x;


    depth = sqrt((radius^2) - (base_line_zerocentered.^2))...
        -radius+height+depth_start;
    
    ref_line5 = [ref_line5, depth];
    
    % Supports
    top_support = zeros(100, 3);
    top_support(:, 1) = top(1);
    top_support(:, 2) = top(2);
    line = linspace(reference_depth, reference_depth+arch_elevation);
    top_support(:, 3) = line';
    
    bottom_support = zeros(100, 3);
    bottom_support(:, 1) = bottom(1);
    bottom_support(:, 2) = bottom(2);
    line = linspace(reference_depth, reference_depth+arch_elevation);
    bottom_support(:, 3) = line';

    
   

    idx_line5 = [15, 23, 24, 25, 26, 36, 41, 42, 43, 44]; 
    coords_line5 = [x_coords(idx_line5), y_coords(idx_line5),...
        z_coords(idx_line5)];

    [~, D] = dsearchn(ref_line5, coords_line5);

    small_arch_dist_RMS = rms(D);

    small_arch_dist_sigma = std(D);
    
    
    % All 
    all_lines = [ref_line1; ref_line2; ref_line3; ref_line4;...
        ref_line5];
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

    plot3(ref_line5(:, 1), ref_line5(:, 2), ref_line5(:, 3),...
        'b', 'LineWidth', 2)

    plot3(top_support(:, 1), top_support(:, 2), top_support(:, 3),...
        'b', 'LineWidth', 2)
    plot3(bottom_support(:, 1), bottom_support(:, 2), bottom_support(:, 3),...
        'b', 'LineWidth', 2)
    axis equal

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
    plot3(ref_line5(:, 1), ref_line5(:, 2), ref_line5(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)

    plot3(top_support(:, 1), top_support(:, 2), top_support(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)
    plot3(bottom_support(:, 1), bottom_support(:, 2), bottom_support(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)

    x_coords = G_updated.Nodes.Coordinates(:, 2)*mm_per_pixel;
    y_coords = 19-G_updated.Nodes.Coordinates(:, 1)*mm_per_pixel;
    z_coords =  0.5*(max(planeNumber)-G_updated.Nodes.Z_coordinates);

    plot(G_updated, 'XData', x_coords, 'YData', y_coords,...
        'ZData', z_coords, 'LineStyle',...
        '-','NodeLabel',{},...
        'EdgeColor', 'b', 'NodeColor', 'r',...
        'LineWidth', 1.5, 'MarkerSize', 6)
    axis equal

    legend('Ground truth','' ,'', '', '', '', '', '3D Reconsruction',...
        'FontSize', 16)

    title("Ground truth model and 3D reconstruction together")

    disp(strcat("------------------ Reconstuction error for prototype #",...
        string(prototype_number), " ----------------"))
    disp(strcat("Filament 1 (Parallell beam on right hand side): RMSE = ", ...
        string(round(parallel_line2_dist_RMS, 2)), " (std = ",...
        string(round(parallel_line2_dist_sigma, 2)), ")"))
    disp(strcat("Filament 2 (Parallell beam on left hand side): RMSE = ", ...
        string(round(parallel_line1_dist_RMS, 2)), " (std = ",...
        string(round(parallel_line1_dist_sigma, 2)), ")"))
    disp(strcat("Filament 3 (Bridge bream): RMSE = ", ...
        string(round(horizontal_line_dist_RMS, 2)), " (std = ",...
        string(round(horizontal_line_dist_sigma, 2)), ")"))
    disp(strcat("Filament 4 (Diagonal bream): RMSE = ", ...
        string(round(diag_beam_dist_RMS, 2)), " (std = ",...
        string(round(diag_beam_dist_sigma, 2)), ")"))
    disp(strcat("Filament 5 (Small arch): RMSE = ", ...
        string(round(small_arch_dist_RMS, 2)), " (std = ",...
        string(round(small_arch_dist_sigma, 2)), ")"))
    disp(strcat("Average: RMSE = ", ...
        string(round(all_dist_RMS, 2)), " (std = ",...
        string(round(all_dist_sigma, 2)), ")"))

elseif prototype_number == 5
    mm_per_pixel = 2/137;

    reference_depth = 4;
    top_elevation = 2; % From mid axis to top is 2.5. Subtracting radius is 2.
    
    % Parallel beam left
    bottom = [352, 848];
    top = [286, 58];
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

    idx_line1 = [2, 3, 4, 5, 6, 7, 8, 21, 22, 23]; % 1 outside
    coords_line1 = [x_coords(idx_line1), y_coords(idx_line1),...
        z_coords(idx_line1)];

    [~, D] = dsearchn(ref_line1, coords_line1);

    parallel_line1_dist_RMS = rms(D);

    parallel_line1_dist_sigma = std(D);

    % Other parallel beam

    bottom = [811, 822];
    top = [759, 34];

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

    idx_line2 = [32, 31, 30, 14, 33, 34, 35, 20, 29, 28, 27]; 
    coords_line2 = [x_coords(idx_line2), y_coords(idx_line2),...
        z_coords(idx_line2)];

    [~, D] = dsearchn(ref_line2, coords_line2);

    parallel_line2_dist_RMS = rms(D);

    parallel_line2_dist_sigma = std(D);

    % Perpendicular bridge
    bottom = [302, 215];
    top = [763, 184];

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

    idx_line3 = [4, 9, 10, 11, 12, 13, 14]; 
    coords_line3 = [x_coords(idx_line3), y_coords(idx_line3),...
        z_coords(idx_line3)];

    [~, D] = dsearchn(ref_line3, coords_line3);

    horizontal_line_dist_RMS = rms(D);

    horizontal_line_dist_sigma = std(D);


    % Elevated bridge
    bottom = [349, 849];
    top = [811, 822];

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

    depth = (reference_depth+top_elevation)*ones(size(ref_line4, 1), 1);

    ref_line4 = [ref_line4, depth];

    idx_line4 = [23, 24, 25, 26, 27]; % Including points with rapid
%     change
%     idx_line4 = [24, 25, 26]; % Only good indices
    coords_line4 = [x_coords(idx_line4), y_coords(idx_line4),...
        z_coords(idx_line4)];

    [~, D] = dsearchn(ref_line4, coords_line4);

    elevated_line_dist_RMS = rms(D);

    elevated_line_dist_sigma = std(D);



    % Supports
    top_support = zeros(100, 3);
    top_support(:, 1) = top(1);
    top_support(:, 2) = top(2);
    line = linspace(reference_depth, reference_depth+top_elevation);
    top_support(:, 3) = line';
    
    bottom_support = zeros(100, 3);
    bottom_support(:, 1) = bottom(1);
    bottom_support(:, 2) = bottom(2);
    line = linspace(reference_depth, reference_depth+top_elevation);
    bottom_support(:, 3) = line';

    % Inclinated bridge
    bottom = [321, 479];
    top = [782, 443];

    bottom = [bottom(1)*mm_per_pixel, 19-bottom(2)*mm_per_pixel];

    top = [top(1)*mm_per_pixel, 19-top(2)*mm_per_pixel];

    point1 = bottom;
    point2 = top;
    n = 100; % Number of points

    % Create n equally spaced values between 0 and 1
    t = linspace(0, 1, n);

    % Preallocate the array for efficiency
    ref_line5 = zeros(n, 2);

    % Calculate the points
    for i = 1:n
        ref_line5(i, :) = point1 + t(i) * (point2 - point1);
    end

    depth = reference_depth + ...
        (hypot(ref_line4(:, 1)-ref_line4(1, 1),...
        ref_line4(:, 2)-ref_line4(1, 2))*(top_elevation/6.7504));

    ref_line5 = [ref_line5, depth];

    idx_line5 = [8, 15, 16, 17, 18, 19, 20]; 
    coords_line5 = [x_coords(idx_line5), y_coords(idx_line5),...
        z_coords(idx_line5)];

    [~, D] = dsearchn(ref_line5, coords_line5);

    inclinated_line_dist_RMS = rms(D);

    inclinated_line_dist_sigma = std(D);

    % Inclinated support
    incl_support = zeros(100, 3);
    incl_support(:, 1) = top(1);
    incl_support(:, 2) = top(2);
    line = linspace(reference_depth, reference_depth+top_elevation);
    incl_support(:, 3) = line';

    
    % All 
    all_lines = [ref_line1; ref_line2; ref_line3; ref_line4;...
        ref_line5];
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
    
    plot3(ref_line5(:, 1), ref_line5(:, 2), ref_line5(:, 3),...
        'b', 'LineWidth', 2)

    % Plot supports
    plot3(incl_support(:, 1), incl_support(:, 2), incl_support(:, 3),...
        'b', 'LineWidth', 2)
    plot3(top_support(:, 1), top_support(:, 2), top_support(:, 3),...
        'b', 'LineWidth', 2)
    plot3(bottom_support(:, 1), bottom_support(:, 2), bottom_support(:, 3),...
        'b', 'LineWidth', 2)
    
    axis equal


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
    plot3(ref_line5(:, 1), ref_line5(:, 2), ref_line5(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)
    
    % Plot supports
    plot3(incl_support(:, 1), incl_support(:, 2), incl_support(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)
    plot3(top_support(:, 1), top_support(:, 2), top_support(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)
    plot3(bottom_support(:, 1), bottom_support(:, 2), bottom_support(:, 3),...
        'Color', '[0.4, 0.4, 0.4]', 'LineWidth', 3)

    x_coords = G_updated.Nodes.Coordinates(:, 2)*mm_per_pixel;
    y_coords = 19-G_updated.Nodes.Coordinates(:, 1)*mm_per_pixel;
    z_coords =  0.5*(max(planeNumber)-G_updated.Nodes.Z_coordinates);

    plot(G_updated, 'XData', x_coords, 'YData', y_coords,...
        'ZData', z_coords, 'LineStyle',...
        '-','NodeLabel',{},...
        'EdgeColor', 'b', 'NodeColor', 'r',...
        'LineWidth', 1.5, 'MarkerSize', 6)
    axis equal

    legend('Ground truth','' ,'', '', '', '', '', '', '3D Reconsruction',...
        'FontSize', 16)

    disp(strcat("------------------ Reconstuction error for prototype #",...
        string(prototype_number), " ----------------"))
    disp(strcat("Filament 1 (Parallell beam on right hand side): RMSE = ", ...
        string(round(parallel_line2_dist_RMS, 2)), " (std = ",...
        string(round(parallel_line2_dist_sigma, 2)), ")"))
    disp(strcat("Filament 2 (Parallell beam on left hand side): RMSE = ", ...
        string(round(parallel_line1_dist_RMS, 2)), " (std = ",...
        string(round(parallel_line1_dist_sigma, 2)), ")"))  
    disp(strcat("Filament 3 (Bridge bream): RMSE = ", ...
        string(round(horizontal_line_dist_RMS, 2)), " (std = ",...
        string(round(horizontal_line_dist_sigma, 2)), ")"))
    disp(strcat("Filament 4 (Diagonal beam): RMSE = ", ...
        string(round(inclinated_line_dist_RMS, 2)), " (std = ",...
        string(round(inclinated_line_dist_sigma, 2)), ")"))
    disp(strcat("Filament 5 (Elevated bream): RMSE = ", ...
        string(round(elevated_line_dist_RMS, 2)), " (std = ",...
        string(round(elevated_line_dist_sigma, 2)), ")"))
    disp(strcat("Average: RMSE = ", ...
        string(round(all_dist_RMS, 2)), " (std = ",...
        string(round(all_dist_sigma, 2)), ")"))

end











