function [sk_structure, ridgeFilt] = mycelium_detection(Igray,...
    polarity, filter_sigma, theta, minimum_branch_length)

[ridgeFilt, theta_max] = ridgeFilter(Igray, filter_sigma,...
    theta, polarity);

BW = binarizeHyphae(ridgeFilt, theta_max);

BW_clean = BW;

% Skeletonize
sk_structure_w_spores = bwskel(BW_clean, 'MinBranchLength',...
    minimum_branch_length);


sk_structure = sk_structure_w_spores;


end