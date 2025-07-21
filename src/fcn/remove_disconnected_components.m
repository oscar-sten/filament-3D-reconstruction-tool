function [sk_structure] = remove_disconnected_components(sk_structure)

[n, m] = size(sk_structure);

frame = zeros(n, m);

% Set the frame = 1
frame(1:5, 1:m) = 1;
frame(n-5:n, 1:m) = 1;
frame(1:n, 1:5) = 1;
frame(1:n, m-5:m) = 1;


disconnected = imfill(~sk_structure, find(frame), 8);


sk_structure = sk_structure & disconnected;


end