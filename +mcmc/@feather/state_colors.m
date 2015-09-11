function a = state_colors(self)
    state_str = self.state_str;

    if ~iscell(state_str)
        state_str = {state_str}; % Wrap one inside the cell
    end
    
    state_s = {'W','R','1','2','3','4'};
    state_color_index = [2 3 4 5 6 6]; % Indexes for cdata

    %states = {'eo','ec','rem','n1','n2','n3','n2s'};

    a = 2*ones(length(state_str),1);
    if isempty(self.state_str{1})
        warning('No colour information, just making them all the same colour');
        return
    end


    for j = 1:length(state_str)
        if length(state_str{j}) == 1
            s = state_str{j};
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
        end
        a(j) = (state_color_index(strcmp(s,state_s)));
    end

