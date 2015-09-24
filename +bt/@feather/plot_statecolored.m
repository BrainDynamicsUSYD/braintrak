function h_out = plot_statecolored(varargin)
	% Plot a single quantity using state_str for the colours
	f = varargin{1};
	varargin = varargin(2:end);
	if ishandle(varargin{1})
		parent_axis = varargin{1};
		varargin = varargin(2:end);
	else
		parent_axis = gca;
	end
	if length(varargin) == 3
		z = varargin{3};
	end
	y = varargin{2};
	x = varargin{1};
    linewidth = 0.5;
    
	if length(x) ~= f.latest
		error('Size mismatch. XYZ data must have same size as the data container')
	end

	[blocks,block_color] = f.state_blocks;
	h_out = [];

	for j = 1:size(blocks,1)
		idx = blocks(j,1):blocks(j,2);
		if length(varargin) == 3
			h_out(end+1) = plot3(x(idx),y(idx),z(idx),'Color',block_color(j,:),'HitTest','off','LineWidth',linewidth,'Parent',parent_axis);
		else
			h_out(end+1) = plot(x(idx),y(idx),'Color',block_color(j,:),'HitTest','off','LineWidth',linewidth,'Parent',parent_axis);
		end
		hold(parent_axis,'on');
	end