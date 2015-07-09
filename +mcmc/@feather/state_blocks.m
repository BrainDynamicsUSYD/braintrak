function [blocks,block_color,state_numbers] = state_blocks(self)
	% For a given state_str
	% return a matrix of the states and their corresponding colours
	% which can be used as indexes for plotting

	state_number = self.state_colors;

	trans = find(diff(state_number))+1;
	trans = [1; trans(:); self.latest];
	blocks = [trans(1:end-1),trans(2:end)];
    cdata = braintrack_utils.state_cdata();
	block_color = cdata(state_number(blocks(:,1)),:);
	state_numbers = state_number(blocks(:,1));
	