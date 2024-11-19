% Create a function to eliminate duplicates while preserving order
function uniqueCellArrays = uniqueOrderPreserved(cellArrays)
    uniqueCellArrays = {};
    seen = containers.Map('KeyType','char', 'ValueType','logical');
    
    for i = 1:numel(cellArrays)
        currentArray = cellArrays{i};
        currentKey = mat2str(currentArray);
        
        if ~isKey(seen, currentKey)
            seen(currentKey) = true;
            uniqueCellArrays{end + 1} = currentArray;
        end
    end
end

% Obtain unique cell arrays while preserving the order
% uniqueCellArrays = uniqueOrderPreserved(allCellArrays);