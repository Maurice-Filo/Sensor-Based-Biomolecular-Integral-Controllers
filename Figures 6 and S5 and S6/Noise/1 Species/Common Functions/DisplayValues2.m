function txt = DisplayValues2(~, event_obj, x1, y1, W1, param_str1, x2, y2, W2, param_str2)
    % Get the selected point position
    pos = get(event_obj, 'Position');

    % Find the index of the selected point in dataset 1
    idx1 = find(x1 == pos(1) & y1 == pos(2), 1);

    if ~isempty(idx1)
        % Point belongs to dataset 1
        associated_values = W1(idx1, :);
        txt = cell(2 + length(associated_values), 1);
        txt{1} = ['x: ', num2str(pos(1))];
        txt{2} = ['y: ', num2str(pos(2))];
        for i = 1:length(associated_values)
            txt{i+2} = [param_str1{i}, ': ', num2str(associated_values(i))];
        end
    else
        % Find the index of the selected point in dataset 2
        idx2 = find(x2 == pos(1) & y2 == pos(2), 1);

        if ~isempty(idx2)
            % Point belongs to dataset 2
            associated_values = W2(idx2, :);
            txt = cell(2 + length(associated_values), 1);
            txt{1} = ['x: ', num2str(pos(1))];
            txt{2} = ['y: ', num2str(pos(2))];
            for i = 1:length(associated_values)
                txt{i+2} = [param_str2{i}, ': ', num2str(associated_values(i))];
            end
        else
            % Point not found in either dataset
            txt = {'Point not found'};
        end
    end
end