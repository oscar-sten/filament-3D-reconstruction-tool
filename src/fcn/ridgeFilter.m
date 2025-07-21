function [ridgeFilt, theta_max] = ridgeFilter(Igray, sigma, theta, polarity)
% Performs steered ridge filtering of the grayscale image Igray.
% Imputs:
% sigma: list of filter scales
% theta: list of angles (deg)
% polarity: bright or dark.

if strcmp(polarity, 'bright')
    Igray2 = imcomplement(Igray);
else
    Igray2 = Igray;
end

for i = 1:length(theta)
    for j = 1:length(sigma)
        ridgeFiltAll(:,:,(((i-1)*length(sigma))+j)) = ...
            steerGaussFilterOrder2(Igray2,...
            theta(i),sigma(j),false);
    end
end

[ridgeFilt, theta_max_index] = max(ridgeFiltAll,[], 3);

theta_max = theta(theta_max_index);

end