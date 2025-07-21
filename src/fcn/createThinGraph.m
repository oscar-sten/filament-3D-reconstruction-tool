function G_thin = createThinGraph(G_original, subGraphPaths, scaling_factor)



G_thin = graph;

for i = 1:numel(subGraphPaths)
    if ~isempty(subGraphPaths{i})
        % Get indices
        indices = subGraphPaths{i};

        % Choose how many indices to keep
        n_keep = round(scaling_factor*sqrt(length(indices)));


        if n_keep < 2
            n_keep = 2;
        end

        if n_keep > numel(indices)
            n_keep = numel(indices);
        end

        % Pick n_keep equidistant elements form indices
        picked_indices = round(linspace(1, length(indices), n_keep));
        picked_elements = indices(picked_indices);

        for j=1:numel(picked_elements)
            if j<numel(picked_elements)
                weight = 0;
                elements_between = indices(picked_indices(j):picked_indices(j+1));
                for k=1:numel(elements_between)-1
                    edge_idx = findedge(G_original, elements_between(k),...
                        elements_between(k+1));
                    weight = weight + G_original.Edges.Weight(edge_idx);
                end
                G_thin = addedge(G_thin, picked_elements(j),...
                    picked_elements(j+1), weight);
            end
        end
    end
end

end