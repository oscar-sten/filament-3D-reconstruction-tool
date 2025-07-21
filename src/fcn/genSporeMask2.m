function [binSporeMask, idSporeMask] = genSporeMask2(centers, radii, rows, columns, skeleton)


[n_circles, ~] = size(centers);

circleImage = false(rows, columns);
circleImageExpanded = zeros(rows, columns);


sporeLabelMargin = 2;
for i=1:n_circles


    center = centers(i, :);
    center = int64(center);
    radius = int64(radii(i));

    margin = round(0.45*radius);
    patch_range = [(center(2)-radius-margin:center(2)+radius+margin)',...
        (center(1)-radius-margin:center(1)+radius+margin)'];
    xlim = rows;
    ylim = columns;
    out_of_bounds_x = (patch_range(:, 1)<1 |...
        patch_range(:, 1)>xlim);
    out_of_bounds_y = (patch_range(:, 2)<1 |...
        patch_range(:, 2)>ylim);
    patch_range(out_of_bounds_x | out_of_bounds_y, :) = [];


    %% TODO Apply Euler based filling. Fill only the ring where the center is
    sk_patch = skeleton(patch_range(:, 1), patch_range(:, 2));

    sk_patch_pruned = imfill(sk_patch, 'holes');
    sk_patch_pruned = imopen(sk_patch_pruned, strel('disk', 1));
    sk_patch_pruned = bwareaopen(sk_patch_pruned, 100); 
    sk_patch_pruned_dil = imdilate(sk_patch_pruned,...
        strel('disk', sporeLabelMargin));

    
    [x, y] = find(sk_patch_pruned);
    bin_mask_ind = sub2ind([rows, columns],...
        x+double(patch_range(1,1)-1), y+double(patch_range(1,2))-1);
    [x, y] = find(sk_patch_pruned_dil);
    id_mask_ind = sub2ind([rows, columns],...
        x+double(patch_range(1,1)-1), y+double(patch_range(1,2)-1));
    circleImage(bin_mask_ind) = true;
    circleImageExpanded(id_mask_ind) = i;


end
binSporeMask = circleImage;
idSporeMask = circleImageExpanded;


end