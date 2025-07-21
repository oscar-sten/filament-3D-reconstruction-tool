function [pixels_scale_bar, mm_scale_bar, position] =...
    read_scale_bar(scale_bar)


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

% Find indeices of scale bar pixels
[scale_bar_x, scale_bar_y] = find(scale_bar);


x_min = min(scale_bar_x);
x_max = max(scale_bar_x);
y_min = min(scale_bar_y);
y_max = max(scale_bar_y);

position = [x_min, x_max, y_min, y_max];

pixels_scale_bar = y_max-y_min;


end