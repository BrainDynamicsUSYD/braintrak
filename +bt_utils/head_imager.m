function p_smoothed = head_imager(varargin)
	% Possible usage
	% head_imager(x,y,P)
	% head_imager(x,y,P,label)
	% head_imager(x,y,f,P) % For an explorer
	% head_imager(x,y,f,P,label)
	% head_imager(x,y,f,P,label,f_secondary,P_secondary)
	smooth_hf = false;

	plot_secondary = false;
	if nargin < 4 || (nargin == 4 && (isempty(varargin{4}) || iscell(varargin{4})))
		explorer_mode = false;
	else
		explorer_mode = true;
	end

	if ~explorer_mode && nargin < 4
		label = [];
	end


	x = varargin{1};
	y = varargin{2};
	if ~explorer_mode
		P = varargin{3};
		if nargin < 4 
			label = [];
		else
			label = varargin{4};
		end
	else
		f = varargin{3};
		P = varargin{4};
		if nargin >= 5
			label = varargin{5};
		else
			label = [];
		end

		if nargin > 5
			plot_secondary = true;
			f_secondary = varargin{6};
			P_secondary = varargin{7};
			assert(size(P,1)==size(P_secondary,1) && size(P,2)==size(P_secondary,2));
		end
	end

	if ~explorer_mode && (size(P,3) > 1 || numel(P)~=numel(x))
		error('Must provide frequencies if frequency-dependent power is specified');
	end



	Lx = 0.5;
	Ly = 0.5;
	x1 = x-Lx/2;
	y1 = y-Ly/2;

	x = x(:);
	y = y(:);

	rmax = bt.data.electrode_positions()/2; % Special mode to set head radius
	valid = sqrt((x1(:).^2+y1(:).^2))<rmax; % These are the only electrodes to include

	[xg,yg] = meshgrid(linspace(0,Lx,40),linspace(0,Ly,40)); % Grid for background colors

	x = x(valid);
	y = y(valid);

	P = reshape_validate(P,valid,explorer_mode);
	
	if plot_secondary
		P_secondary = reshape_validate(P_secondary,valid,explorer_mode);
	end

	if explorer_mode && smooth_hf
		f_filt = f > 25; % Set this to 5 to smooth the alpha peak and ensure that the frequencies don't change
		
		% THIS IS WORKING CODE FOR INDEPENDENT SMOOTHING
		% for j = 1:size(P,2)
		% 	P(f_filt,j) =smooth(P(f_filt,j),10);
		% end
		% if plot_secondary
		% 	f_filt = f_secondary > 25;
		% 	for j = 1:size(P,2)
		% 		P_secondary(f_filt,j) =smooth(P_secondary(f_filt,j),10);
		% 	end
		% end

		% 3D GRID INTERPOLATION - SMOOTHING IN SPACE AND FREQUENCY
		pp = P(f_filt,:);
		xp = repmat(x(:)',size(pp,1),1);
		yp = repmat(y(:)',size(pp,1),1);
		fp = repmat(f(f_filt),1,size(pp,2));

		fl = linspace(fp(1,1),fp(end,1),50).';
		fl = f(f_filt);
		[xp1,yp1,fp1] = meshgrid(xg(1,:),yg(:,1),fl);
		%pp1 = griddata(xp(:),yp(:),fp(:),pp(:),xp1,yp1,fp1);
		pp_spec = scatteredInterpolant(xp(:),yp(:),fp(:),pp(:));
		pp1 = pp_spec(xp1,yp1,fp1);

		pp1a = smooth3(pp1,'box',[5 5 5]);

		% DIAGNOSTIC PLOT SHOWING SMOOTHED SURFACE
		% figure
		% idx = 20;
		% f_test = fp(idx,1);
		% [~,idx2] = min(abs(fl-f_test));
		% f_test
		% fl(idx2)
		% surf(xp1(:,:,20),yp1(:,:,20),pp1a(:,:,20))
		% hold on
		% P2 = P(f_filt,:);
		% scatter3(x,y,P2(20,:),'ro')

		pp2 = griddata(xp1(:),yp1(:),fp1(:),pp1a(:),xp,yp,fp);
		P(f_filt,:) = pp2;

	end
	p_smoothed = P;

	figure

	if explorer_mode % Do the initialization and all that jazz
		subplot(1,2,1);
	end
	img = imagesc([min(xg(:)) max(xg(:))],[min(yg(:)) max(yg(:))],zeros(size(xg)));
	hold on
	utils.draw_head(rmax,0.4,[Lx/2 Ly/2])
	draw_electrodes(x,y,label)
	axis equal
	set(gca,'XLim',[0 0.5],'YLim',[0 0.5])
	colorbar('SouthOutside')
	img.HitTest = 'off';
	title('Spatial power distribution')
	xlabel('X (m)');
	ylabel('Y (m)');
	hold off

	if explorer_mode % Do the initialization and all that jazz
		ax_1 = gca;
		hold on
		elec_marker = scatter(x(1),y(1),'ro','hittest','off');
		hold off
	
		subplot(1,2,2);
		ax_2 = gca;
		spec = loglog(f,squeeze(P(:,1)));
		spec.HitTest = 'off';
		ax_2.XLim = [1 45];
		ax_2.YLim = [min(P(:)) max(P(:))];
		fb = ax_2.XBaseline;
		fb.Color = 'r';
		fb.BaseValue = 0;
		fb.Visible = 'on';
		title('Power spectrum')
		xlabel('Frequency (Hz)')
		ylabel('Power spectral density')
		if plot_secondary
			hold on
			spec_secondary = loglog(f_secondary,squeeze(P_secondary(:,1)),'r--');
			spec_secondary.HitTest = 'off';
		end
		set(ax_1,'ButtonDownFcn', @(a,b)  update_spec(b) )
		set(ax_2,'ButtonDownFcn', @(a,b)  update_map(b) )
	else
		[Xi,Yi,Zi] = griddata(x,y,P,xg,yg,'cubic'); 
		img.CData = Zi;
	end



	function update_spec(b)
		click_value = b.IntersectionPoint(1:2);
		[~,idx] = min((sum(bsxfun(@minus,[x(:) y(:)],click_value).^2,2)));

		spec.YData = P(:,idx); % ... so they are reversed here into (row,col) = (y,x)
		if plot_secondary
			spec_secondary.YData = P_secondary(:,idx);
		end
		elec_marker.XData = x(idx);
		elec_marker.YData = y(idx);
	end

	function update_map(b)
		f_click = b.IntersectionPoint(1);
		
		[~,f_idx] = min(abs(f-f_click));
		fb.BaseValue = f(f_idx);
		cdata = P(f_idx,:);
		xt = x;
		yt = y;
		if length(unique(x))==1
			xt = [xt(:)-0.05; xt(:)+0.05];
			yt = [yt(:); yt(:)];
			cdata = [cdata(:); cdata(:)];
		end

		if length(unique(yt)) == 1
			yt = [yt(:)-0.05; yt(:)+0.05];
			xt = [xt(:); xt(:)];
			cdata = [cdata(:); cdata(:)];		
		end

		[~,~,Zi] = griddata(xt,yt,cdata(:),xg,yg,'cubic'); 
		img.CData = Zi;
		img.AlphaData = isfinite(Zi);
		ax_1.CLim = [min(cdata(:)) max(cdata(:))];
	end

end

function P2 = reshape_validate(P,valid,explorer_mode)
	if explorer_mode % 2D output
		if length(size(P)) == 2 % The first dimension is frequency
			P2 = zeros(size(P,1),sum(valid));
		else
			P2 = zeros(size(P,3),sum(valid));
		end
	else % 1D output
		P2 = zeros(sum(valid),1);
	end

	count = 1;
	for j = 1:length(valid)
		if valid(j)
			if explorer_mode && length(size(P)) == 2 % 2D output, and it is frequency dependent
				P2(:,count) = P(:,j);
			elseif explorer_mode % 3D output, remap array indexes
				[a,b] = ind2sub([size(P,1) size(P,2)],j);
				P2(:,count) = P(a,b,:);
			else % 1D output
				P2(count) = P(j);
			end
			count = count + 1;
		end
	end
end

function draw_electrodes(x,y,label)

	scatter(x,y,'HitTest','off')

	if ~isempty(label)
		for j = 1:length(label)
			text(x(j),y(j)-0.007,label{j},'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'hittest','off')
		end
	end
end


