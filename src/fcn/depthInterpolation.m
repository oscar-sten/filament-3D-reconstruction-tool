function G_updated = depthInterpolation(G_updated,...
    reliability_threshold, considered_span, imageMetaData)

FM_lists = G_updated.Nodes.FM_list;
planeNumberLists = G_updated.Nodes.FM_list_planes;

planeNumber = imageMetaData.planeNumber;

z_estimates = zeros(numel(FM_lists), 1);

for i=1:numel(FM_lists)
    current_list = FM_lists{i};
    current_plane_numbers = planeNumberLists{i};
    fm_values = current_list;
    x = 1:numel(current_plane_numbers);
    % Remove NaN
    x(isnan(fm_values)) = [];
    fm_values = fm_values(~isnan(fm_values));
    if numel(fm_values)>3
        % A-priori guess
        [~, idx] = max(fm_values);
        idx_interval = idx-considered_span:idx+considered_span;
        idx_interval(idx_interval<1) = [];
        idx_interval(idx_interval>numel(x)) = [];
        estimation_interval_fm = fm_values(idx_interval);
        estimation_interval_x = x(idx_interval);
        try
            [fitObject, ~] = fit(estimation_interval_x',...
                estimation_interval_fm', "gauss1");
            x_range = linspace(min(x), max(x));

            val = feval(fitObject, x_range);
            [~, idx_max] = max(val);
            mu_final = x_range(idx_max);
            A = fitObject.a1;
            s = fitObject.c1;

            G = A.*exp(-((estimation_interval_x-mu_final)./s).^2);

            erf = estimation_interval_fm-G;
            fmax = max(current_list);
            P=numel(erf);
            R = 20*log10(P*fmax/(sum(abs(erf))));
            % R is in dB. R greater than 0 -> more signal than noise.
            % R is analogus to Peak Signal to Noise Ratio (PSNR)
            if R<reliability_threshold
                mu_final = -1;
            end

        catch
            mu_final = -1;
        end

    else
        mu_final = -1;
    end
    z_estimates(i) = mu_final;

end


% Mean value

idx_unknown = find(z_estimates==-1);

disp(strcat("Nodes with insufficient reliability: ", ...
    string(round(((numel(idx_unknown)/numel(z_estimates))*100), 1)), "%"));

G_updated.Nodes.Z_coordinates = z_estimates;



end

