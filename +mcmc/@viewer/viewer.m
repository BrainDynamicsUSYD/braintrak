classdef viewer < handle
	% A viewer should be initialized with all of the data that is known
	properties
		parent_figure

		spec_ax
			spec_exp
			spec_fit
			spec_fit_2
			electrode_titles

		hist_ax
			hist_line
			hist_val
		chisq_ax
			chisq_ax_txt
			chisq_vert_line

		tent_ax
			xyz_marker
			vol3d_data
			red_trail 
			title_str 
			tent_mesh
			selected_electrode_index = 0

		plot_chisq
		plot_xyz_posterior = [false false false]
		n_hist % The total number of histograms - n_params + n_xyz_posteriors
		model
		is_spatial = false
		current_index = 1 % Time index of the data currently on the screen
		feather
		feather_2
	end

	methods

		function self = viewer(varargin)
			% USAGES
			% viewer(feather) - Make an animation
			% viewer(feather,feather_2)
			% viewer.update_struct(fit_data,plot_data)


			self.model = varargin{1}.model;
			self.feather = varargin{1};
			if isempty(self.feather.time)
				self.feather.time = 1:length(self.feather.fit_data);
			end

			if self.feather.latest == 1
				self.plot_chisq = false;
			else
				self.plot_chisq = true;
			end

			if nargin == 2
				self.feather_2 = varargin{2};
				if self.feather_2.latest ~= self.feather.latest
					error('Data objects must have the same number of data')
				end
				if all(size(self.feather_2.fit_data(1).target_P)~=size(self.feather.fit_data(1).target_P))
					error('Fit output size mismatch');
				end
				if all(self.feather_2.fit_data(1).target_P(:,1)~=self.feather.fit_data(1).target_P(:,1))
					warning('Target power spectrum is different for the two data sets - is this intentional?')
				end
			end

			self.parent_figure = figure('Visible','off');

			if isfield(self.feather.fit_data,'xyz_posterior')
				self.plot_xyz_posterior = ~self.model.xyz_present();
				% self.plot_xyz_posterior(1) = ~any(strcmp(self.model.param_names,'X'));
				% self.plot_xyz_posterior(2) = ~any(strcmp(self.model.param_names,'Y'));
				% self.plot_xyz_posterior(3) = ~any(strcmp(self.model.param_names,'Z'));
			end

			% First, work out how many rows and columns are required
			n_cols = 4;
			n_param_rows = ceil((self.model.n_params+sum(self.plot_xyz_posterior))/2); % Total number of params to plot
			if ~isprop(self.model,'electrodes')
				% Old fits have only one electrode
				self.model.set_electrodes('Main')
			end

			if isa(self.model,'mcmc.model.template_spatial')
				self.is_spatial = true;
			end

			n_spec_rows = ceil(length(self.model.electrodes)/2);
			% Give the tent the top two rows on the LHS
			% Then, the bottom panel requires n_spec_rows
			lhs_rows = 2+n_spec_rows;
			rhs_rows = n_param_rows; 
			n_rows = max(lhs_rows,rhs_rows);
			a_idx = @(a,b) n_cols*(a-1) + b; % This function gives back the subplot number for row/col a/b

			% Proposed grid is n_cols x n_rows
			% Now, decide how much to allocate to each
			if length(self.model.electrodes) == 1
				% Divide them up more or less equally
				tent_row_count = ceil(n_rows/2);
				spec_row_count = floor(n_rows/2);
			else
				% First allocate the 
				lhs_excess = lhs_rows - n_spec_rows;
				tent_row_count = 2; % Could do something smarter here, but not going to
				spec_rows_count = n_spec_rows;
			end

			% Initialize the axes
			tent_ax_pos = [];
			for j = 1:tent_row_count
				for k = 1:(n_cols - 2)
					tent_ax_pos(end+1) = a_idx(j,k);
				end
			end

			spec_ax_pos = [];
			if length(self.model.electrodes) == 1
				for j = 1+tent_row_count:n_rows
					for k = 1:(n_cols - 2)
						spec_ax_pos(end+1) = a_idx(j,k);
					end
				end
				self.spec_ax = subplot(n_rows+self.plot_chisq,n_cols,spec_ax_pos,'Parent',self.parent_figure);
				title(self.model.electrodes{1})
			else
				count = 1;
				for j = 1+tent_row_count:n_rows
					for k = 1:(n_cols - 2)
						if count <= length(self.model.electrodes)
							self.spec_ax(end+1) = subplot(n_rows+self.plot_chisq,n_cols,a_idx(j,k),'Parent',self.parent_figure);
							self.electrode_titles(end+1) = title(self.model.electrodes{count});
							count = count + 1;
						end
					end
				end
			end

			for j = (n_rows/2+1):n_rows
				for k = 1:(n_cols - 2)
					spec_ax_pos(end+1) = a_idx(j,k);
				end
			end

			% ---------------------------
			self.tent_ax = subplot(n_rows+self.plot_chisq,n_cols,tent_ax_pos,'Parent',self.parent_figure);
			xyz = [0 0 0];
			hold(self.tent_ax,'on');
			self.xyz_marker = scatter3(self.tent_ax,xyz(1),xyz(2),xyz(3),50,'go','MarkerFaceColor','g');
			self.vol3d_data = [];
			self.red_trail = plot3(self.tent_ax,xyz(1),xyz(2),xyz(3),'r','LineWidth',0.5);
			self.title_str = title('TENT');
			self.tent_mesh = tent.surface_hg2;
			self.tent_mesh.EdgeColor = 'flat';
			self.tent_mesh.Visible = 'off';

			% DRAWING COLORBAR DISABLES ROTATION IN R2014b
			%cbar = colorbar(self.tent_ax);
			%ylabel(cbar,'Frequency (Hz)');

			self.tent_ax.CLim = [0 60];
			self.tent_ax.ALim = [0 1];
			hold(self.tent_ax,'off');
			axis(self.tent_ax,'equal');
			axis(self.tent_ax,'vis3d');
			set(self.tent_ax,'DataAspectRatio',[0.75 1 1]);
			self.tent_ax.XTick = 0:0.5:1;
			self.tent_ax.YTick = -1:0.5:1;
			self.tent_ax.ZTick = 0:0.5:1;
			box(self.tent_ax,'on')
			grid(self.tent_ax,'on')
			xlabel('X');
			ylabel('Y');
			zlabel('Z');


			% ----------------------------
			% Initialize spectra
			
			for j = 1:length(self.model.electrodes)
				set(self.spec_ax(j),'XScale','log','YScale','log')
				hold(self.spec_ax(j),'on');
				box(self.spec_ax(j),'on');
				set(self.spec_ax(j),'XLim',[0.1 45]);
				%set(self.spec_ax(j),'YLim',[3e-4 1]);
				xtf = [1:2:10 15 20 30 45];
				set(self.spec_ax(j),'XTick',xtf);

				self.spec_exp(j) = plot(self.spec_ax(j),1,1,'HitTest','off');
				self.spec_fit(j) = plot(self.spec_ax(j),1,1,'r','HitTest','off');
				if ~isempty(self.feather_2)
					self.spec_fit_2(j) = plot(self.spec_ax(j),1,1,'--','Color',[0 0.5 0],'HitTest','off');
				end
				hold(self.spec_ax(j),'off');
				%self.spec_ax(j).UserData = j; % The userdata stores the electrode index
				if self.is_spatial
					set(self.spec_ax(j),'ButtonDownFcn',@(a,b,c) self.set_spec_electrode(j))
				end
			end

			% ----------------------------

			% Initialize histograms
			self.hist_ax = [];
			self.n_hist = self.model.n_params + sum(self.plot_xyz_posterior);
			xyz_strs = {'X','Y','Z'};
			xyz_units = {'','',''};

			param_symbols = [self.model.param_symbols,xyz_strs(self.plot_xyz_posterior)]; % Append only the relevant title strs
			param_units = [self.model.param_units,xyz_units(self.plot_xyz_posterior)]; % Append only the relevant title strs

			lim = self.model.limits;
			xyz_lim = [0  -1     0 ; 1   1   1.3]; % It is a bit ugly that these appear in multiple places
			lim = [lim xyz_lim(:,self.plot_xyz_posterior)];

			for j = 1:self.n_hist
				self.hist_ax(end+1) = subplot(n_rows+self.plot_chisq,n_cols,a_idx(ceil(j/2),3+mod(j-1,2)),'Parent',self.parent_figure);
				box(self.hist_ax(j),'on');
				hold(self.hist_ax(j),'on');
				set(self.hist_ax(j),'YLim',[0 1]);
				xlim_vals = [lim(1,j) lim(2,j)];
				xlim_vals(xlim_vals == eps | xlim_vals == -eps) = 0;
				set(self.hist_ax(j),'XLim',xlim_vals);
				self.hist_line(j) = plot(self.hist_ax(j),[1 1],[1 1],'b');
				self.hist_val(j) = plot(self.hist_ax(j),1,1,'r','LineWidth',2);
				hold(self.hist_ax(j),'off');
				if isempty(param_units{j})
					title(self.hist_ax(j),param_symbols{j})
				else
					title(self.hist_ax(j),sprintf('%s (%s)',param_symbols{j},param_units{j}))
				end
			end

			% Initialize chisq
			if self.plot_chisq
				chisq_ax_pos = [];
				for j = 1:n_cols
					chisq_ax_pos(end+1) = a_idx(n_rows+1,j);
				end
				self.chisq_ax = subplot(n_rows+1,n_cols,chisq_ax_pos,'Parent',self.parent_figure);
				box(self.chisq_ax,'on')
				hold(self.chisq_ax,'on')
				self.chisq_ax_txt = utils.fig_letter('1');
				self.chisq_ax_txt.FontSize = 14;
				self.chisq_vert_line = plot(self.chisq_ax,[1 1],[0 0],'k','HitTest','off');
				hold(self.chisq_ax,'off')
				self.prep_animation()
			else
				self.update_idx(1);
			end

			self.parent_figure.Visible = 'on';
		end

		function update_struct(self,fit_data,plot_data,title_str,xyz,fit_data_2)
			if nargin < 6 || isempty(fit_data_2)
				fit_data_2 = [];
			end

			if nargin < 5 || isempty(xyz)
				xyz = fit_data.xyz;
			end

			if nargin < 4 || isempty(title_str)
				if ~isfield(fit_data,'state_str')
					title_str = '';
				else
					title_str = fit_data.state_str;
				end
			end

			% Update the tent
			set(self.tent_ax,'XLim',plot_data.xyzlim(:,1),'YLim',plot_data.xyzlim(:,2),'ZLim',[0 1]);

			p = self.model.p_from_params(fit_data.fitted_params);

			if self.is_spatial && self.selected_electrode_index ~= 0
				p = p.params_at_xy(self.model.output_x(self.selected_electrode_index),self.model.output_y(self.selected_electrode_index));
				set(self.xyz_marker,'XData',p.xyz(1),'YData',p.xyz(2),'ZData',p.xyz(3));
				fit_data.fitted_params = self.model.params_from_p(p);
			else
				set(self.xyz_marker,'XData',fit_data.xyz(end,1),'YData',fit_data.xyz(end,2),'ZData',fit_data.xyz(end,3));
				set(self.red_trail,'XData',xyz(:,1),'YData',xyz(:,2),'ZData',xyz(:,3));
			end

			[x,y,z,u] = tent.compute(p.alpha(1),p.beta(1),p.t0,p.gammae);
			edge_alphadata = x+y<0.90 & x > 0 & x < 1;
			set(self.tent_mesh,'XData',x,'YData',y,'ZData',z,'CData',u,'AlphaData',+edge_alphadata,'Visible','on');
			self.vol3d_data = self.vol3d_alpha_fast(self.tent_ax,plot_data.xyzlim,+(plot_data.Vcount>0),plot_data.Vcount,self.vol3d_data);
			self.title_str.String = title_str;

			% Update the spectrum
			for j = 1:size(fit_data.target_P,2)
				set(self.spec_exp(j),'XData',fit_data.target_f,'YData',fit_data.target_P(:,j));
				set(self.spec_fit(j),'XData',fit_data.target_f,'YData',fit_data.fitted_P(:,j));
				if ~isempty(fit_data_2)
					set(self.spec_fit_2(j),'XData',fit_data_2.target_f,'YData',fit_data_2.fitted_P(:,j));
				end
			end

			% Update the histograms
			idx1 = find(self.plot_xyz_posterior);
			for j = 1:self.n_hist
				if j <= self.model.n_params
					%normconst = max(fit_data.posterior_pp.y(:,j));
					normconst = 1;
					set(self.hist_line(j),'XData',fit_data.posterior_pp.x(:,j),'YData',fit_data.posterior_pp.y(:,j)/normconst);
					ymax = max(get(self.hist_line(j),'YData'));
					set(self.hist_val(j),'XData',[fit_data.fitted_params(j) fit_data.fitted_params(j)],'YData',[0 ymax]);
					set(self.hist_ax(j),'YLim',[0 ymax]);

					if fit_data.skip_fit(j)
						set(self.hist_line(j),'LineStyle','--');
					else
						set(self.hist_line(j),'LineStyle','-');
					end

					if fit_data.skip_fit(j) == 2
						set(self.hist_val(j),'LineStyle','--');
					else
						set(self.hist_val(j),'LineStyle','-');
					end
				else
					idx = idx1(j-self.model.n_params);
					normconst = 1;
					set(self.hist_line(j),'XData',fit_data.xyz_posterior.x(:,idx),'YData',fit_data.xyz_posterior.y(:,idx)/normconst);
					ymax = max(get(self.hist_line(j),'YData'));
					set(self.hist_val(j),'XData',[fit_data.xyz(idx) fit_data.xyz(idx)],'YData',[0 ymax]);
					set(self.hist_ax(j),'YLim',[0 ymax]);
				end
			end

			drawnow update
		end

		function set_spec_electrode(self,idx)

			if self.selected_electrode_index ~= 0
				set(self.electrode_titles(self.selected_electrode_index),'Color','k');
			end

			if self.selected_electrode_index == idx
				self.selected_electrode_index = 0;
			else
				set(self.electrode_titles(idx),'Color','r');
				self.selected_electrode_index = idx;
			end

			[a,b] = self.feather.retrieve(self.current_index);
			if ~isempty(self.feather_2)
				a2 = self.feather_2.retrieve(idx);
			else
				a2 = [];
			end
			xyz = self.feather.xyz;
			self.update_struct(a,b)
			%self.update_struct(a,b,[],xyz(max(1,self.current_index-100):self.current_index,:),a2);
		end

		function prep_animation(self)
			chisq = self.feather.chisq;
			contaminated = braintrack_utils.chisq_outliers(chisq);
			chisq(contaminated) = NaN;
			
			fixed_idx = arrayfun(@(x) x.skip_fit(1),self.feather.fit_data) == 3;
			time_idx = self.feather.time;

			axes(self.chisq_ax);
			hold(self.chisq_ax,'on');
			self.feather.plot_statecolored(time_idx,chisq);
			scatter(time_idx(fixed_idx),chisq(fixed_idx),'r.');
			hold(self.chisq_ax,'off');
			axis(self.chisq_ax,'tight');
			set(self.chisq_ax,'XLim',[time_idx(1) time_idx(end)]);
			yl = get(self.chisq_ax,'YLim');
			self.chisq_vert_line.YData = [yl(1) yl(2)];

			fcn = @(idx) self.update_idx(idx);
			
			set(self.parent_figure,'WindowButtonMotionFcn',@(d,e,f) fcn([]))
			set(self.chisq_ax,'ButtonDownFcn',@(a,b,c) playback_controller(self.parent_figure,fcn,self.feather.latest));

			plot_persistent.idx = 1;
			plot_persistent.play = 0;
			plot_persistent.mouseover = false;
			plot_persistent.vert_line = self.chisq_vert_line;
			guidata(self.parent_figure,plot_persistent);

			fcn(1);
		end

		function update_idx(self,idx)
			% Update all the plots to a specific index
			if utils.is_multiple_call()
				return;
			end

			pdata = guidata(self.parent_figure);
			if ~isempty(pdata) && ~pdata.mouseover && isempty(idx)
				return
			end

			if isempty(idx)
				t = get(self.chisq_ax,'CurrentPoint');
				% Now find the time closest to it - since this axis corresponds to time
				[~,idx] = min(abs(self.feather.time-t(1)));

				% idx = round(idx(1));
				% if idx < 1
				% 	idx = 1;
				% elseif idx > self.feather.latest
				% 	idx = self.feather.latest;
				% end

				pdata.idx = idx;
				guidata(self.parent_figure,pdata);
			end

			[a,b] = self.feather.retrieve(idx);
			if ~isempty(self.feather_2)
				a2 = self.feather_2.retrieve(idx);
			else
				a2 = [];
			end
			xyz = self.feather.xyz;
			if self.plot_chisq
				self.chisq_vert_line.XData = [self.feather.time(idx) self.feather.time(idx)];
				self.chisq_ax_txt.String = num2str(self.feather.time(idx));
			end
			self.update_struct(a,b,[],xyz(max(1,idx-100):idx,:),a2);
			self.current_index = idx;
		end	
	end

	methods (Static)

		function cache = vol3d_alpha_fast(parent,xyzlim,cdata,alphadata,cache)
		    % Volume render 3-D data with transparency 
		    % Feed the model back in to re-render the plot
		    % vol3d_alpha(p)
		    % vol3d_alpha(limits,cdata,alphadata)
		    if nargin < 5 || isempty(cache)
		        cache.h_out = zeros(1,sum(size(cdata)));
		        cache.h_present = zeros(1,sum(size(cdata)));
		        cache.h_zeroed = ones(1,sum(size(cdata)));
		    end

		    if nargin < 4 || isempty(alphadata)
		        alphadata = ones(size(cdata));
		    else
		        alphadata = alphadata/max(alphadata(:));
		    end

		    if isempty(xyzlim)
		        xyzlim = [0 0 0; 1 1 1];
		    end

		  %  cdata = +(cdata>0);

		    gridres = diff(xyzlim(:,3))/size(cdata,3);
		    xdata = xyzlim(1,1):gridres:xyzlim(2,1);
		    ydata = xyzlim(1,2):gridres:xyzlim(2,2);
		    zdata = xyzlim(1,3):gridres:xyzlim(2,3);

		    h_idx = 0;
		    for dim = 1:3 % For each dimension
		        for j = 1:size(cdata,dim) % For each slice of the dimension
		            h_idx = h_idx + 1;
		            [cdat,adat] = mcmc.viewer.fetch_slice(cdata,alphadata,dim,j);

		            if ~cache.h_present(h_idx) && any(isfinite(adat(:)))
		                b = mcmc.viewer.fetch_xyz(xdata,ydata,zdata,dim,j);
		                cache.h_out(h_idx) = surface(b{1},b{2},b{3},'Parent',parent,'facecolor','texture','edgealpha',0,'facealpha','texturemap'); % Drawing the surface using only the endpoints gives a huge performance boost
		                cache.h_present(h_idx) = 1;
		            end

		            if any(isfinite(adat(:)))
		                adat(~isfinite(cdat)) = 0;
		                set(cache.h_out(h_idx),'CData',cdat,'alphadata',adat);
		                cache.h_zeroed(h_idx) = 0;
		            elseif cache.h_present(h_idx) && ~cache.h_zeroed(h_idx)
		                set(cache.h_out(h_idx),'alphadata',zeros(size(cdat)));
		                cache.h_zeroed(h_idx) = 1;
		            end
		        end
		    end
		   
		end

		function a = fetch_xyz(xdata,ydata,zdata,dim,j);
		    switch dim
		        case 1
		            a{1} = xdata([1,end;1,end]);
		            a{2} = ydata(j*ones(2));
		            a{3} = zdata([1,1;end,end]);
		        case 2
		            a{1} = xdata(j*ones(2));
		            a{2} = ydata([1,end;1,end]);
		            a{3} = zdata([1,1;end,end]);
		        case 3
		            a{1} = xdata([1,end;1,end]);
		            a{2} = ydata([1,1;end,end]);
		            a{3} = zdata(j*ones(2));
		    end    
		end

		function [cdat,adat] = fetch_slice(cdata,alphadata,dim,j);
		    switch dim
		        case 1
		            cdat = reshape(cdata(j,:,:),size(cdata,2),size(cdata,3))';
		            adat = reshape(alphadata(j,:,:),size(alphadata,2),size(alphadata,3))';
		        case 2
		            cdat = reshape(cdata(:,j,:),size(cdata,1),size(cdata,3))';
		            adat = reshape(alphadata(:,j,:),size(alphadata,1),size(alphadata,3))';
		        case 3
		            cdat = reshape(cdata(:,:,j),size(cdata,1),size(cdata,2));
		            adat = reshape(alphadata(:,:,j),size(alphadata,1),size(alphadata,2));
		    end    
		end
	end
end



% If user left clicks, playback starts and mouseover disabled
% If user right clicks, mouseover is toggled. 
function playback_controller(h,fcn,stop_idx)
	if strcmp(get(gcf,'SelectionType'),'open')
		return
	end

	switch get(gcf,'SelectionType')
		case 'open'
			return
		case 'alt'
			p = guidata(h);
			if ~p.mouseover
				p.mouseover = true;
				set(p.vert_line,'Color','r')
				p.play = 0;
			else
				set(p.vert_line,'Color','k');
				p.mouseover = false;
			end
			guidata(h,p);
		case 'normal'
			p = guidata(h);
			if p.play
				p.play = 0;
				p.mouseover = false;
				set(p.vert_line,'Color','k');
				guidata(h,p);
				return;
			else
				p.play = 1;
				p.mouseover = false;
				set(p.vert_line,'Color','k');

				guidata(h,p);
				for j = p.idx:stop_idx
					if ~mod(j-1,10)
						drawnow
					end
					p = guidata(h);
					if p.play
						fcn(j);
					else
						p.play = 0; % Guarantee playback is stopped
						p.idx = j;
						guidata(h,p);
						return
					end
				end
			end
		
	end
end