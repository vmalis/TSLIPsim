function closestIndex = findClosestElementIndex(array, targetValue)
    % Calculate the absolute differences between targetValue and each element
    absoluteDifferences = abs(array - targetValue);
    
    % Find the index of the element with the minimum absolute difference
    [~, closestIndex] = min(absoluteDifferences);
end
