function h_out = scatter_statecolored(varargin)
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
    
	if length(x) ~= f.latest
		error('Size mismatch. XYZ data must have same size as the data container')
	end

	cdata = bt_utils.state_cdata;
	cdata = cdata(f.state_colors,:);

	if length(varargin) == 3
		h_out = scatter3(x,y,z,10,cdata,'Marker','.','HitTest','off','Parent',parent_axis);
	else
		h_out = scatter(x,y,10,cdata,'Marker','.','HitTest','off','Parent',parent_axis);
	end
