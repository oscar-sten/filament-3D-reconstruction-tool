% Shallow crossing detection in 2D

clear all
rng('default')
addpath(genpath("fcn\"));
addpath(genpath("..\data\"))

% Sample selection
sample = "s.jpg"; % Available choices include: "s.jpg", "o.jpg", "b.jpg",
% "p.jpg", "r.jpg", and "v.jpg".

%% Parameters

filter_sigma = 4;
spore_diameter_range_mm = [0.04 0.165];
theta = [0:15:360];
polarity = 'dark';
patch_range = [10 10];
minimum_branch_length = 20;
densityFactor = 1;
side_margin = 10;

overlap_displacement_tolerance = 0; % Displacement tolerance in number
% of pixels in each direction.

%% Read image

current_image = imread(sample);

current_image = rgb2gray(current_image);

[current_image, pixels_scale_bar, mm_scale_bar] = remove_scale_bar_2(...
        current_image);

mm_per_pixel = mm_scale_bar/pixels_scale_bar;

Img = imread(sample);

%% Process skeleton structure

[sk_structure, ridgeFilt] = mycelium_detection(current_image,...
    polarity, filter_sigma, theta, minimum_branch_length);

%% Remove disconnected components
[sk_structure] = remove_disconnected_components(sk_structure);


%% Transform into representative/thin graph
[G_thin, G_minimal] = get_G_thin(sk_structure, densityFactor);


%% In 2D no merging is done, all nodes derive from the same plane.
G_thin.Nodes.Z_coordinates = ones(numel(G_thin.Nodes.Z_planes), 1);

G_thin.Nodes.SourceGraphIdx = ones(numel(G_thin.Nodes.Z_planes), 1);
%% Identify overlaps

G_input = G_thin;

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
plot(G_input, 'XData', x_coords, 'YData', y_coords, 'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b')

n_test = size(crossings, 1);

for i=1:2:n_test-1
    crossings_tmp = crossings(i:i+1, :);
    plot(crossings_tmp(:, 2), crossings_tmp(:, 1), 'LineWidth',...
        3, 'Color', color_list(mod(i, numel(color_list))+1))
end
title("Identified overlaps")


% Coordinates of perfect crossings, useful for evaluation
degrees = degree(G_minimal);
nodes_degree_4 = find(degrees == 4); % nodes of degree 4
nodes_degree_4IDs = G_minimal.Nodes.Name(nodes_degree_4);
perfect_crossing_coords = G_minimal.Nodes.Coordinates(nodes_degree_4, :);

%% Match overlaps
G_updated = match_overlaps(G_input, G_minimal, nodeIDCrossingPairs);

% Visualize
x_coords = G_updated.Nodes.Coordinates(:, 2);
y_coords = G_updated.Nodes.Coordinates(:, 1);

figure,
imshow(Img)
hold on 
plot(G_updated,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b')
title("Graph with overlaps")


%% Compute morphological graph with overlaps

% This is useful for evaluation
degree3Nodes = find(degree(G_updated) == 3);
BP_nodes = G_updated.Nodes.Name(degree3Nodes);


% Minimal morphological graph with overlaps
G_updated2 = thin2MorphGraph(G_updated);



x_coords = G_updated2.Nodes.Coordinates(:, 2);
y_coords = G_updated2.Nodes.Coordinates(:, 1);

figure,
imshow(Img)
hold on 
plot(G_updated2,'XData', x_coords,'YData',y_coords,'LineStyle',...
    '-','NodeLabel',{},...
    'EdgeColor', 'r', 'NodeColor', 'b');


%% Validation test

idx_BP = findnode(G_updated, BP_nodes);
BP_coords = G_updated.Nodes.Coordinates(idx_BP, :);

n_cross = size(nodeIDCrossingPairs, 1);

crossing_centers = zeros(n_cross, 2);

for i=1:n_cross
    idx = findnode(G_thin, string(nodeIDCrossingPairs(i, :)'));
    coords = G_thin.Nodes.Coordinates(idx, :);
    center = mean(coords, 1);
    crossing_centers(i, :) = center;
end

crossing_centers = [crossing_centers; perfect_crossing_coords];


figure,
imshow(Img)
hold on 
plot(BP_coords(:, 2), BP_coords(:, 1), '^', 'LineWidth', 2, 'Color', ...
    "[0.4940 0.1840 0.5560]")
plot(crossing_centers(:, 2), crossing_centers(:, 1), 's',...
    'LineWidth', 2, 'Color', "b")

legend('Identified Branch Points', 'Identified Overlap Points')


%% Compare to GT

test = convertStringsToChars(sample);
suffix_annotation = '_intersections_label_special.mat';

Img = imread(sample);

GT = load(strcat(test(1), suffix_annotation));
intersection_labels = GT.gTruth.LabelData;


figure,
imshow(Img)
hold on

% Define sets
% Branches
all_branches = [BP_coords(:, 2), BP_coords(:, 1)];
identified_branches = all_branches;
correctly_identified_branches = [];
false_negative_branches = [];
% Crosssings
identified_crossings = [crossing_centers(:, 2), crossing_centers(:, 1)];
correctly_identified_crossings = [];
false_negative_crossings = [];


% Mycelium detection errors
detection_errors = table2array(intersection_labels(:, "Detection_Error"));
detection_errors = detection_errors{1};
n_detection_errors = size(detection_errors, 1);

for i=1:n_detection_errors
    currentBox = detection_errors(i, :);
    x_interval = [currentBox(1), currentBox(1)+currentBox(3)];
    y_interval = [currentBox(2), currentBox(2)+currentBox(4)];

    idx = find((identified_branches(:, 1)>x_interval(1)) & ...
        identified_branches(:, 1)<x_interval(2) & ...
        identified_branches(:, 2)>y_interval(1) & ...
        identified_branches(:, 2)<y_interval(2));
    if ~isempty(idx)
        identified_branches(idx, :) = [];
    end

    idx = find((identified_crossings(:, 1)>x_interval(1)) & ...
        identified_crossings(:, 1)<x_interval(2) & ...
        identified_crossings(:, 2)>y_interval(1) & ...
        identified_crossings(:, 2)<y_interval(2));
    if ~isempty(idx)
        identified_crossings(idx, :) = [];
    end

    rectangle('Position', detection_errors(i, :), 'EdgeColor','k');
end

% Quantify accuracy

% Rectangle coordinates are in fromat [x y w h] where x and y define the
% upper left corner and w is width and h is height.

% Branches

branchBoxes = table2array(intersection_labels(:, 'Branch'));
branchBoxes = branchBoxes{1};

n_branches = size(branchBoxes, 1);

for i=1:n_branches
    currentBox = branchBoxes(i, :);
    x_interval = [currentBox(1), currentBox(1)+currentBox(3)];
    y_interval = [currentBox(2), currentBox(2)+currentBox(4)];
    idx = find((identified_branches(:, 1)>x_interval(1)) & ...
        identified_branches(:, 1)<x_interval(2) & ...
        identified_branches(:, 2)>y_interval(1) & ...
        identified_branches(:, 2)<y_interval(2));
    if ~isempty(idx)
        correct_branch = identified_branches(idx, :);
        correctly_identified_branches = [correctly_identified_branches;...
            correct_branch];
        identified_branches(idx, :) = [];
    else
        false_negative_branches = [false_negative_branches;...
            [mean(x_interval), mean(y_interval)]];
    end
    rectangle('Position', branchBoxes(i, :), 'EdgeColor','[0, 0.6, 0]');
end

plot(correctly_identified_branches(:, 1),...
    correctly_identified_branches(:, 2), 'o', 'Color', 'k')

if ~isempty(false_negative_branches)
    plot(false_negative_branches(:, 1), false_negative_branches(:, 2), '^', ...
        'Color', '[0.9290 0.6940 0.1250]', 'LineWidth', 2)
end

if ~isempty(identified_branches)
    plot(identified_branches(:, 1), identified_branches(:, 2), '^', ...
        'Color', 'red', 'LineWidth', 2)
end


% Crossings/Overlaps

crossingBoxes = table2array(intersection_labels(:, "Crossing"));
crossingBoxes = crossingBoxes{1};
n_crossings = size(crossingBoxes, 1);

for i=1:n_crossings
    currentBox = crossingBoxes(i, :);
    x_interval = [currentBox(1), currentBox(1)+currentBox(3)];
    y_interval = [currentBox(2), currentBox(2)+currentBox(4)];
    idx = find((identified_crossings(:, 1)>x_interval(1)) & ...
        identified_crossings(:, 1)<x_interval(2) & ...
        identified_crossings(:, 2)>y_interval(1) & ...
        identified_crossings(:, 2)<y_interval(2));
    if ~isempty(idx)
        correct_crossing = identified_crossings(idx, :);
        correctly_identified_crossings = [correctly_identified_crossings;...
            correct_crossing];
        identified_crossings(idx, :) = [];
    else
        false_negative_crossings = [false_negative_crossings;...
            [mean(x_interval), mean(y_interval)]];
    end
    rectangle('Position', crossingBoxes(i, :), 'EdgeColor','b');
end


plot(correctly_identified_crossings(:, 1),...
    correctly_identified_crossings(:, 2), 'o', 'Color', 'k')

if ~isempty(false_negative_crossings)
    plot(false_negative_crossings(:, 1), false_negative_crossings(:, 2),...
        's', 'Color', '[0.9290 0.6940 0.1250]', 'LineWidth', 2)
end

if ~isempty(identified_crossings)
    plot(identified_crossings(:, 1), identified_crossings(:, 2), ...
        's', 'Color', 'red', 'LineWidth', 2)
end


n_tp_branch = size(correctly_identified_branches, 1);
n_fp_branch = size(identified_branches, 1);
n_fn_branch = size(false_negative_branches, 1);

n_tp_crossing = size(correctly_identified_crossings, 1);
n_fp_crossing = size(identified_crossings, 1);
n_fn_crossing = size(false_negative_crossings, 1);

figure()
groups = categorical(["Branches", "Crossings"]);
b = bar(groups, [n_tp_branch, n_fp_branch, n_fn_branch;...
    n_tp_crossing, n_fp_crossing, n_fn_crossing], 'grouped');

b(1).FaceColor = [0, 0.6, 0];
b(2).FaceColor = 'red';
b(3).FaceColor = [0.9290 0.6940 0.1250];

% Put text above bars
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips3 = b(3).XEndPoints;
ytips3 = b(3).YEndPoints;
labels3 = string(b(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')


legend('True positive', 'False positive', 'False negative')

%% Precision
precision_branches = n_tp_branch/(n_tp_branch + n_fp_branch);
precision_crossings = n_tp_crossing/(n_tp_crossing + n_fp_crossing);

% Recall
recall_branches = n_tp_branch/(n_tp_branch + n_fn_branch);
recall_crossings = n_tp_crossing/(n_tp_crossing + n_fn_crossing);

% F1 score
f1_branches = (2*n_tp_branch)/(2*n_tp_branch + n_fn_branch + n_fp_branch);
f1_crossings = (2*n_tp_crossing)/(2*n_tp_crossing +...
    n_fn_crossing + n_fp_crossing);

% Print stats

disp(strcat("---------------------- Stats of test image: ", sample,...
    " ----------------------"))
disp("----------------------- Branches --------------------------")

disp(strcat("True positives, branches: ", string(n_tp_branch)))

disp(strcat("False positives, branches: ", string(n_fp_branch)))

disp(strcat("False negatives, branches: ", string(n_fn_branch)))

disp(strcat("Precision in identifying branches: ", ...
    string(round(precision_branches*100, 2)), "%"))

disp(strcat("Recall in identifying branches: ", ...
    string(round(recall_branches*100, 2)), "%"))

disp(strcat("F1 score on identifying branches: ", ...
    string(round(f1_branches, 2))))

disp("----------------------- Overlaps --------------------------")

disp(strcat("True positives, branches: ", string(n_tp_crossing)))

disp(strcat("False positives, branches: ", string(n_fp_crossing)))

disp(strcat("False negatives, branches: ", string(n_fn_crossing)))

disp(strcat("Precision in identifying overlaps: ", ...
    string(round(precision_crossings*100, 2)), "%"))

disp(strcat("Recall in identifying overlaps: ", ...
    string(round(recall_crossings*100, 2)), "%"))

disp(strcat("F1 score on identifying overlaps: ", ...
    string(round(f1_crossings, 2))))



