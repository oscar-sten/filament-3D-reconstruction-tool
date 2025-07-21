function [patch, x_range_1, y_range_1] = extractPatch(I, coords, patch_range)


% Image size.
[n, m] = size(I);

x = coords(1);
y = coords(2);

% Define patch indices
x_range = x-patch_range(1):x+patch_range(2);
y_range = y-patch_range(1):y+patch_range(2);

% Remove invalid indices
x_range(x_range<1)=[];
x_range(x_range>n)=[];
y_range(y_range<1)=[];
y_range(y_range>m)=[];

% Extract patch
patch = I(x_range, y_range);

x_range_1 = x_range(1);
y_range_1 = y_range(1);
end