function [image, pixels_scale_bar, mm_scale_bar, threshold] = remove_scale_bar_2(I)
% Author Oscar Sten, IIT, 10/11/2021
% Fully automatic
% Identify the scale bar, which is typically either completely black or
% completely white

% Assuming the scalebar is black and that the unit is either mm or um.
% TODO: Make sure that the whole scale bar is included
threshold = 1;

while 1
    scale_bar = (I<threshold);
    
    % Read scale bar with OCR
    OCR = ocr(scale_bar);
    
    % See the number of words
    n_words = numel(OCR.Words);
    if n_words > 2
        % error('Too many words identified')
        number_read = 0;
        unit_read = 0;
    elseif n_words < 1
        %         error('No text identified')
        number_read = 0;
        unit_read = 0;
    elseif n_words == 1
        word = OCR.Words{1};
        for i=1:numel(word)
            if isempty(str2num(word(i)))
                break;
            end
        end
        number = word(1:i-1);
        unit = word(i:numel(word));
    else
        number = OCR.Words{1};
        unit = OCR.Words{2};
    end
    
    % Identify the number
    try
        if (str2num(number(1)) == 0) && (str2num(number)>1)
            % It's common to miss the point, but in this case it can be solved
            % quite easily.
            value = str2num(number);
            n_digits=numel(num2str(value));
            value = value*(10^-n_digits);
        else
            value = str2num(number);
        end
        number_read = 1;
    catch
        number_read = 0;
    end
    
    if number_read == 1
        % Identify the unit and return value in mm
        try
            if unit == 'mm'
                mm_scale_bar = value;
            elseif ((unit(1) == 'p')||(unit(1) == 'u')) && (unit(2) == 'm')
                mm_scale_bar = value*0.001;
            else
                error(strcat('Invalid unit: ', unit))
            end
            unit_read = 1;
        catch
            unit_read = 0;
        end
    end
    
    if (number_read == 1) && (unit_read == 1)&& ...
            ((sum(OCR.WordConfidences>0.75)/...
            length(OCR.WordConfidences))==1)
        break;
    else
        threshold = threshold+1;
    end
    
    if threshold > 255
        error('The text could not be identified with the desired confidence')
    end
    
end
 
    
% Find indeices of scale bar pixels
[scale_bar_x, scale_bar_y] = find(scale_bar);
x_min = min(scale_bar_x);
x_max = max(scale_bar_x);
y_min = min(scale_bar_y);
y_max = max(scale_bar_y);

% Take the pixels from the are representing a 10 pixels wide frame around
% the scale bar area
above = I(x_min:x_max, (y_min-20):(y_min-10));


% Calculate the mean color from the frame area and replace the scalebar
% with the mean surrounding color

mean_color = mean(above(:));
I((x_min-10):(x_max+10), (y_min-10):(y_max+10)) = mean_color;

image = I;
pixels_scale_bar = y_max-y_min;
end