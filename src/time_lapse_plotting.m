% This script reproduces the graphs in Figure 10.


clear all
rng('default')

addpath(genpath("fcn\"));


data_path_prefix1 = "..\data\Imgs4TimeLapse\Imgs4TimeLapse1\";
data_path_prefix2 = "..\data\Imgs4TimeLapse\Imgs4TimeLapse2\";

data_suffix = ".tif";


% Filaments do not result with the same polarity in the different planes.
% Filaments in focus result bright, and slightly out of focus result dark,
% this is because the time lapse was acquired with a Nikon ECLIPSE Ti2 
% inverted microscope in brightfield mode. For this reason, the filament
% detection and the overlap detection steps are adapted, and only one focal
% plane is considered.

% Parameters
filter_sigma = 3;
D_T = 12;
patch_range = [10 10];
minimum_branch_length = 5;
theta = 0:15:360;
lower_lim = 6;
sensitivity = 0.8;
method = "PhaseCode";

structuring_element = strel('disk', 3);


mm_per_pixel = 1.48*0.001;

spore_diameter_range_mm = [0.025 0.14];

lower_lim_pix = (spore_diameter_range_mm(1)/2)/mm_per_pixel;
upper_lim_pix = (spore_diameter_range_mm(2)/2)/mm_per_pixel;

scale_factor = lower_lim_pix/lower_lim;

uppper_lim = ceil(upper_lim_pix/scale_factor);

patch_size_mm = 0.05;
patch_size = round(patch_size_mm/mm_per_pixel);


z_plane_number = "19";

polarity = "dark";

name_base = "time_lapse_hyphae";


time_list = ["01", "61"];

data_path_prefix = data_path_prefix1;


total_length_list = zeros(1, 26);
total_length_list2 = zeros(1, 26);
N_edges_list = zeros(1, 26);
N_tips_list = zeros(1, 26);
N_bp_list = zeros(1, 26);

for i=1:2

    time_instance = time_list(i);

    image_name = strcat(name_base, "t", time_instance, "z", z_plane_number);

    img_path = strcat(data_path_prefix, image_name, data_suffix);

    Img = imread(img_path);

    Igray = rgb2gray(Img);

    Igray = Igray(450:1342, 450:1300);

    [rows, columns] = size(Igray);


    [ridgeFilt2, theta_max] = ridgeFilter(Igray, filter_sigma,...
        theta, polarity);

    BW = binarizeHyphae(ridgeFilt2, theta_max);

    BW2 = imdilate(BW, structuring_element);

    BW3 = bwmorph(BW2, 'majority');

    BW3_skel = bwskel(BW3, 'MinBranchLength', 5);

    sk_structure = BW3_skel;


    current_image_rescaled = imresize(Igray, 1/scale_factor);
    [centers, radii, metric] = imfindcircles(current_image_rescaled,...
        [lower_lim uppper_lim], "ObjectPolarity", "dark",...
        "Sensitivity", sensitivity, "Method", method);

    centers = scale_factor*centers;
    radii = scale_factor*radii;

    [binSporeMask, idSporeMask] = genSporeMask2(centers, radii, rows,...
        columns, sk_structure);

    sk_structure = bwareaopen(sk_structure, 200);
    sk_structure = sk_structure & ~binSporeMask;
    sk_structure = bwareaopen(sk_structure, 50);


    G_thin = get_G_thin(sk_structure, 1);

    G_thin = simplify(G_thin); 

    G_updated = identify_overlaps_special(Igray, G_thin,...
        filter_sigma, strel('disk', 2));

    G_updated = simplify(G_updated);

    G_updated2 = thin2MorphGraph(G_updated);

    % Extract relevant data
    total_length_list(i) = mm_per_pixel*sum(G_updated2.Edges.Weight); 
    N_edges_list(i) = numel(G_updated2.Edges.Weight);
    N_tips_list(i) = numel(find(degree(G_updated2)==1));
    N_bp_list(i) = numel(find(degree(G_updated2)==3));

end


%% djskjsk

time_list = ["007", "013", "019", "025", "031", "037", "043", "048",...
    "055", "061", "066", "073", "079", "085", "091", "097", "103", ...
    "109", "115", "121", "127", "133", "139", "145"];

data_path_prefix = data_path_prefix2;

for i=1:24

    time_instance = time_list(i);

    image_name = strcat(name_base, "t", time_instance, "z", z_plane_number);

    img_path = strcat(data_path_prefix, image_name, data_suffix);

    Img = imread(img_path);

    Igray = rgb2gray(Img);

    Igray = Igray(450:1342, 450:1300);

    [rows, columns] = size(Igray);
    

    [ridgeFilt2, theta_max] = ridgeFilter(Igray, filter_sigma,...
        theta, polarity);

    BW = binarizeHyphae(ridgeFilt2, theta_max);

    BW2 = imdilate(BW, structuring_element);

    BW3 = bwmorph(BW2, 'majority');


    BW3_skel = bwskel(BW3, 'MinBranchLength', 5);

    sk_structure = BW3_skel;


    current_image_rescaled = imresize(Igray, 1/scale_factor);
    [centers, radii, metric] = imfindcircles(current_image_rescaled,...
        [lower_lim uppper_lim], "ObjectPolarity", "dark",...
        "Sensitivity", sensitivity, "Method", method);

    centers = scale_factor*centers;
    radii = scale_factor*radii;

    [binSporeMask, idSporeMask] = genSporeMask2(centers, radii, rows,...
        columns, sk_structure);

    sk_structure = bwareaopen(sk_structure, 200);
    sk_structure = sk_structure & ~binSporeMask;
    sk_structure = bwareaopen(sk_structure, 50);



    G_thin = get_G_thin(sk_structure, 1);

    G_thin = simplify(G_thin); 

    G_updated = identify_overlaps_special(Igray, G_thin,...
        filter_sigma, strel('disk', 2));

    G_updated = simplify(G_updated);

    G_updated2 = thin2MorphGraph(G_updated);


    % Extract relevant data
    total_length_list(i+2) = mm_per_pixel*sum(G_updated2.Edges.Weight);
    N_edges_list(i+2) = numel(G_updated2.Edges.Weight);
    N_tips_list(i+2) = numel(find(degree(G_updated2)==1));
    N_bp_list(i+2) = numel(find(degree(G_updated2)==3));

end




%% Actual plotting


figure,
plot(0:25, total_length_list, '-*', 'LineWidth', 2, 'Color', 'b')
xlabel("Hours", "FontSize", 14, 'Interpreter', 'latex')
ylabel("Total length [mm]", "FontSize", 14, 'Interpreter', 'latex')


figure,
plot(0:25, N_edges_list, '-*', 'LineWidth', 2, 'Color', 'r')
xlabel("Hours", "FontSize", 14, 'Interpreter', 'latex')
ylabel("Number of filaments", "FontSize", 14, 'Interpreter', 'latex')



figure,
plot(0:25, N_tips_list, '-*', 'LineWidth', 2, 'Color', 'k')
xlabel("Hours", "FontSize", 14, 'Interpreter', 'latex')
ylabel("Number of tips", "FontSize", 14, 'Interpreter', 'latex')


figure,
plot(0:25, N_bp_list, '-*', 'LineWidth', 2, 'Color', '[0.8500 0.3250 0.0980]')
xlabel("Hours", "FontSize", 14, 'Interpreter', 'latex')
ylabel("Number of branches", "FontSize", 14, 'Interpreter', 'latex')



