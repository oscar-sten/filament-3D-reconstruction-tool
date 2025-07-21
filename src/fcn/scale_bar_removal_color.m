function [I, mm_per_pixel, position] = scale_bar_removal_color(I)


scale_bar = (I(:, :, 2) > I(:, :, 1));

[pixels_scale_bar, mm_scale_bar, position] = read_scale_bar(scale_bar);

I = rgb2gray(I);

above = I(position(1)-60:position(2)-60,...
    (position(3)):(position(4)));


mean_color = mean([above(:)]);

I((position(1)-10):(position(2)+10),...
    (position(3)-10):(position(4)+10)) = mean_color;

mm_per_pixel = mm_scale_bar/pixels_scale_bar;
end