function traversals = traversalsFromPairs(nodeIDCrossingPairs)
% Input: A set of filament ID pairs that are connected through a crossing
% Output: Cell array with filament ID traversals.
% Algorithm 1 in manuscript of https://doi.org/10.1016/j.ecoinf.2024.102670

traversals = {};

if ~isempty(nodeIDCrossingPairs)
    U = nodeIDCrossingPairs(:, 1);
    V = nodeIDCrossingPairs(:, 2);
else
    U = [];
    V = [];
end

while ~isempty(U)
    traversal = [];
    currentID = U(1);
    traversal = [traversal; currentID];
    while ismember(currentID, U) || ismember(currentID, V)
        X = find(currentID == U);
        Y = find(currentID == V);
        if ~isempty(X)
            correspID = V(X(1));
            U(X(1)) = [];
            V(X(1)) = [];
            traversal = [traversal; correspID];

        elseif ~isempty(Y)
            correspID = U(Y(1));
            U(Y(1)) = [];
            V(Y(1)) = [];
            traversal = [traversal; correspID];

        end
    currentID = correspID;
    end

    traversals{end+1} = traversal;

end

end