function [width, color] = width_measure(patch, method)

% Input: patch : [N x N], method : 'Otsu'
% Output: filament width in pixels

if strcmp(method, 'Otsu')
    bin_patch = ~imbinarize(mat2gray(patch));
else
    error('Invalid method')
end

% Apply morphological majority filter
bin_patch = bwmorph(bin_patch, 'majority');

% Measure filament width along it's medial axis
eucl_dist = bwdist(~bin_patch);

bin_patch_skel = bwskel(bin_patch);

width_values = 2*eucl_dist(bin_patch_skel);

width = mean(width_values);

color = mean(patch(bin_patch));

end