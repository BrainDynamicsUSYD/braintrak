function cones(f)
	    %if ~isnumeric(varargin{1})
	    	%f = varargin{1};
    	xyz = f.xyz;
    	chisq = f.chisq;
    	contaminated = braintrack_utils.chisq_outliers(chisq);
    	fitted = isfinite(chisq) & ~contaminated;
    	valid = [fitted(2:end)==1 & fitted(1:end-1)==1]; % This says 
	    % else
	    %     xyz = varargin{1};
	    %     xyzdiff = varargin{2};
	    %     valid = true(size(xyz,1),1);
	    % end
	    
	    
	    xyzdiff = diff(xyz); 
	    xyz = xyz(1:end-1,:);
        save temp2 xyz xyzdiff

        % load temp
        % valid = ones(size(xyz,1),1);

        xyzlim = [0  -1     0 ;...
    	          1   1   1.3];

        gridres = 0.05;


        xlim = (xyz(:,1) < xyzlim(2,1) & xyz(:,1) > xyzlim(1,1));
        ylim = (xyz(:,2) < xyzlim(2,2) & xyz(:,2) > xyzlim(1,2));
        zlim = (xyz(:,3) < xyzlim(2,3) & xyz(:,3) > xyzlim(1,3));

        xyz = xyz(xlim(:) & ylim(:) & zlim(:) & valid(:),:); 
        xyzdiff = xyzdiff(xlim(:) & ylim(:) & zlim(:) & valid(:),:);


        xBinIndex = floor((xyz(:,1)-xyzlim(1,1))./gridres)+1;
        yBinIndex = floor((xyz(:,2)-xyzlim(1,2))./gridres)+1;
        zBinIndex = floor((xyz(:,3)-xyzlim(1,3))./gridres)+1;


        % Arrays for the total velocity
        xv = zeros((xyzlim(2,2)-xyzlim(1,2))/gridres,(xyzlim(2,1)-xyzlim(1,1))/gridres,(xyzlim(2,3)-xyzlim(1,3))/gridres);
        yv = xv;
        zv = xv;
        Vcount = xv;

        for zval = 1:max(zBinIndex)
    	    this_z = zBinIndex==zval;
    	
    	    xb = xBinIndex(this_z);
    	    yb = yBinIndex(this_z);
    	    this_diff = xyzdiff(this_z,:);

    	    for j = 1:length(xb)-1
    	    	xv(yb(j),xb(j),zval) =  xv(yb(j),xb(j),zval)+ this_diff(j,1);
    	    	yv(yb(j),xb(j),zval) =  yv(yb(j),xb(j),zval)+ this_diff(j,2);
    	    	zv(yb(j),xb(j),zval) =  zv(yb(j),xb(j),zval)+ this_diff(j,3);

    		    Vcount(yb(j),xb(j),zval) =  Vcount(yb(j),xb(j),zval)+1;
    	    end
        end

        xv = xv./Vcount;
        yv = yv./Vcount;
        zv = zv./Vcount;

        xvals = (xyzlim(1,1)+gridres:gridres:xyzlim(2,1))-gridres/2;
        yvals = (xyzlim(1,2)+gridres:gridres:xyzlim(2,2))-gridres/2;
        zvals = (xyzlim(1,3)+gridres:gridres:xyzlim(2,3))-gridres/2;

		[xg,yg,zg] = meshgrid(xvals,yvals,zvals);

		% filt = Vcount(:)>0;
		% xv = xv(filt);
		% yv = yv(filt);
		% zv = zv(filt);
		% xg = xg(filt);
		% yg = yg(filt);
		% zg = zg(filt);


		figure
		%scatter3(xg1(~points_present),yg1(~points_present),zg1(~points_present),1,'b.')
		%scatter3(xg1(points_present),yg1(points_present),zg1(points_present),50,'g.')
		%quiver3(xg,yg,zg,xdir,ydir,zdir,2);

		hcones = coneplot(xg,yg,zg,xv,yv,zv,xg,yg,zg,0,'nearest')
		hold on
		set(hcones,'FaceColor','red','EdgeColor','none')
		se.paper_surface_hg2

function [xyz_out,xyzdiff] = load_all
    xyzdiff = [];
    xyz_out = [];
    for j = 1:9
        f = mcmc.feather.import(sprintf('postmaster/control_apnea_art_%d_postmaster',j));

        xyz = f.xyz;
        chisq = f.chisq;
        contaminated = braintrack_utils.chisq_outliers(chisq);
        fitted = isfinite(chisq) & ~contaminated;
        valid = [fitted(2:end)==1 & fitted(1:end-1)==1]; % This says 

        xyzdiff = [xyzdiff; diff(xyz)]; 
        xyz_out = [xyz_out; xyz(1:end-1,:)];
    end