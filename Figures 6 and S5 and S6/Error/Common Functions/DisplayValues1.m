function txt = DisplayValues1(~, event_obj, x1, y1, W1, param_str1)
    % Get the selected point position
    pos = get(event_obj, 'Position');
    
    % Find the index of the selected point in dataset 1
    idx1 = find(x1 == pos(1) & y1 == pos(2), 1);
    
    % Point belongs to dataset 1
    associated_values = W1(idx1, :);
    txt = cell(2 + length(associated_values), 1);
    txt{1} = ['x: ', num2str(pos(1))];
    txt{2} = ['y: ', num2str(pos(2))];
    for i = 1:length(associated_values)
        txt{i+2} = [param_str1{i}, ': ', num2str(associated_values(i))];
    end
end