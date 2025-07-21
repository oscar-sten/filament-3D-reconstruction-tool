function G_updated = smoothDepthCoords(G_updated)

z_estimates = G_updated.Nodes.Z_coordinates;
idx_unknown = find(z_estimates==-1);
for i=1:numel(idx_unknown)
    idx = idx_unknown(i);
    node_id = G_updated.Nodes.Name(idx);
    nghbrs = neighbors(G_updated, node_id);
    nghbrZCoords = G_updated.Nodes.Z_coordinates(...
        findnode(G_updated, nghbrs));
    nghbrZCoords(nghbrZCoords==-1) = [];

    if numel(nghbrZCoords)>0
        est_depth = mean(nghbrZCoords);
    else
        % If the immediate neighbors are nullvalued, 
        % find the closests neighbor that isn't.
        explored = [];
        while numel(nghbrZCoords)<1
            nghbrs = setdiff(nghbrs, explored);
            explored = [explored; nghbrs];
            nghbrs_tmp = [];
            for j=1:numel(nghbrs)
                nghbrs2 = neighbors(G_updated, nghbrs{j});
                nghbrs2 = setdiff(nghbrs2, explored);
                nghbrZCoords = G_updated.Nodes.Z_coordinates(...
                    findnode(G_updated, nghbrs2));
                nghbrZCoords(nghbrZCoords==-1) = [];
                nghbrs_tmp = [nghbrs_tmp; nghbrs2];
            end
            nghbrs = nghbrs_tmp;
            
            if isempty(nghbrs)
                % If there are no valid neighbors available we simply
                % consider the MAP-option.
                nghbrZCoords = G_updated.Nodes.SourceGraphIdx;
                break;
            end
        end
        est_depth = mean(nghbrZCoords);
    end
    G_updated.Nodes.Z_coordinates(idx) = est_depth;
end