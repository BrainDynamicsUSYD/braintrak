function clouds(f)
	figure
	
	xyz = f.xyz;
	state_color = f.state_colors;
	cdata = braintrak_utils.state_cdata;



    xyzlim = [0  -1     0 ;...
	          1   1   1.3];

    gridres = 0.025;

    for j = 1:size(cdata,1)
		pts = state_color == j;
		if any(pts)
			[Vcount] = utils.simple_bin3d(xyz(pts,:),xyzlim,gridres);
			vol = mcmc.viewer.vol3d_alpha_fast(gca,xyzlim,j*+(Vcount>0),Vcount);
			hold on
		end
	end
	colormap(cdata)
	set(gca,'CLim',[1,size(cdata,1)])
	tent.surface_hg2