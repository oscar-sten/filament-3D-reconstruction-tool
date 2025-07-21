function lmax = nonMaxSuppression(magMax, theta_max)
% Perform Non Maximum Suppression on the filter response matrix magMax in
% the corresonding directions in theta_max
% Inputs:
% magMax [n x m], filter responses on interest.
% theta_max [n x m], corresponding angles (deg) yielding the maximal filter
% response.
% Output:
% lmax [n x m], binary matrix where 1 indicates that a pixel ins a local
% maxima in the direction yielding the maximal ridge filter response (i.e.,
% the direction perpendicular to the ridge).

% Adapted from:
% Rachmawan (2022). Canny Edge Detection 
% (https://www.mathworks.com/matlabcentral/fileexchange/46859-canny-edge-detection), 
% MATLAB Central File Exchange. Retrieved March 17, 2022


[n, m] = size(magMax);

main_dir=zeros(n, m);

%Adjusting directions to nearest 0, 45, 90, or 135 degree
for i = 1  : n
    for j = 1 : m
        if ((theta_max(i, j) >= 0 ) && (theta_max(i, j) < 22.5) || ...
                (theta_max(i, j) >= 157.5) && (theta_max(i, j) < 202.5) ||...
                (theta_max(i, j) >= 337.5) && (theta_max(i, j) <= 360))
            main_dir(i, j) = 0;
        elseif ((theta_max(i, j) >= 22.5) && (theta_max(i, j) < 67.5) || ...
                (theta_max(i, j) >= 202.5) && (theta_max(i, j) < 247.5))
            main_dir(i, j) = 45;
        elseif ((theta_max(i, j) >= 67.5 && theta_max(i, j) < 112.5) ||...
                (theta_max(i, j) >= 247.5 && theta_max(i, j) < 292.5))
            main_dir(i, j) = 90;
        elseif ((theta_max(i, j) >= 112.5 && theta_max(i, j) < 157.5) ||...
                (theta_max(i, j) >= 292.5 && theta_max(i, j) < 337.5))
            main_dir(i, j) = 135;
        end
    end
end

BW = zeros(n,m);

%Non-Maximum Supression
for i=2:n-1
    for j=2:m-1
        if (main_dir(i,j)==0)
            BW(i,j) = (magMax(i,j) == max([magMax(i,j),...
                magMax(i,j+1), magMax(i,j-1)]));
        elseif (main_dir(i,j)==45)
%             BW(i,j) = (magMax(i,j) == max([magMax(i,j),...
%                 magMax(i+1,j-1), magMax(i-1,j+1)]));
            BW(i,j) = (magMax(i,j) == max([magMax(i,j),...
                magMax(i+1,j+1), magMax(i-1,j-1)]));

        elseif (main_dir(i,j)==90)
            BW(i,j) = (magMax(i,j) == max([magMax(i,j),...
                magMax(i+1,j), magMax(i-1,j)]));
        elseif (main_dir(i,j)==135)
%             BW(i,j) = (magMax(i,j) == max([magMax(i,j),...
%                 magMax(i+1,j+1), magMax(i-1,j-1)]));
              BW(i,j) = (magMax(i,j) == max([magMax(i,j),...
                magMax(i+1,j-1), magMax(i-1,j+1)]));
        end
    end
end

lmax = BW;

end