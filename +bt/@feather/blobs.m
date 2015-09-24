function blobs(f)
	figure
	
	xyz = f.xyz;
	state_color = f.state_colors;
	cdata = braintrak_utils.state_cdata;


	for j = 1:size(cdata,1)
		pts = state_color == j;
		if any(pts)
			h = utils.alphavol_wrapper(xyz(pts,:),Inf);

			hold on
			set(h,'FaceColor',cdata(j,:),'EdgeColor','none','FaceAlpha',0.9,'UserData','blob');
		end
	end

	se.paper_surface_hg2