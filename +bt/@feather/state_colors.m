function a = state_colors(self)
    % Convert the state string into a colour

    state_str = self.state_str;
    [~,state_names,short_mapping] = bt_utils.state_cdata;

    if ~iscell(state_str)
        state_str = {state_str}; % Wrap one inside the cell
    end

    
    a = ones(length(state_str),1);

    for j = 1:length(state_str)
        if ~isempty(state_str{j})
            if length(state_str{j}) == 1
                a(j) = short_mapping{2}(find(strcmp(state_str{j},short_mapping{1})));
            else
                idx = find(strcmp(lower(state_str{j}),lower(state_names)));
                if ~isempty(idx)
                    a(j) = idx;
                else
                    d = sscanf(state_str{j},'%*s (%d/%d)');
                    if isempty(d)
                        error('Empty! Did you pass the state string AND run compress()?')
                    end

                    if d(1) >= d(2)
                        s = state_str{j}(1);
                    else
                        s = state_str{j}(3);
                    end
                    a(j) = short_mapping{2}(find(strcmp(s,short_mapping{1})));
                end
            end
        end
    end

