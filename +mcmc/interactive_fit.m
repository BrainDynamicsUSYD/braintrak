function interactive_fit(fit)
	h.m = fit.model;
	h.pars_orig = fit.fitted_params;
	h.f = fit.fit_data.target_f;
	h.f = linspace(1,45,1000);
	h.fig = figure;


	h.ax = axes('Parent',h.fig,'Position',[0.07 0.45 0.45 0.5]);
	[f,P] = h.m.spectrum(h.pars_orig,h.f);
	s1 = loglog(f,P);
	xlabel('Frequency (Hz)');
	ylabel('Power (normalized)')
	hold on
	h.s2 = loglog(f,P,'r--');
	set(h.ax,'XLim',[1 45]);
	loglog(fit.fit_data.target_f,fit.fit_data.target_P,'b--')

	h.t0_marker(1) = plot([NaN NaN],[NaN NaN],'g--')
	h.t0_marker(2) = plot([NaN NaN],[NaN NaN],'g--')
	h.t0_marker(3) = plot([NaN NaN],[NaN NaN],'g--')
	h.spindle_marker = plot([NaN NaN],[NaN NaN],'m--')
	h.spindle_text = text(NaN,NaN,'\sigma','VerticalAlignment','bottom')
	h.t0_text(1) = text(NaN,NaN,'f_\alpha','VerticalAlignment','bottom')
	h.t0_text(2) = text(NaN,NaN,'2f_\alpha','VerticalAlignment','bottom')
	h.t0_text(3) = text(NaN,NaN,'3f_\alpha','VerticalAlignment','bottom')
	hold off

	h.tent_ax = axes('Parent',h.fig,'Position',[0.6 0.35 0.35 0.6]);
	se.paper_surface_hg2;
	delete(findobj(h.tent_ax,'type','Surface'))
	[x,y,z,u] = tent.compute;
	h.tent_surf = mesh(x,y,z,u);
	h.tent_marker = scatter3(NaN,NaN,NaN,50,'go','MarkerFaceColor','g');
	h.tent_drop = plot3([NaN NaN],[NaN NaN],[NaN NaN],'r--');
	colorbar('Location','southoutside')
	set(h.tent_ax,'CLim',[0 45])
	n = h.m.n_params;
	col = 2;
	row = 2*ceil(n/2); % Slider row and label

	rowheight = 0.3/row; % Normalized units per row
	row_offset = 0:rowheight:(row-1)*(rowheight);

	set(gcf,'Color',[0.94 0.94 0.94])
	for j = 1:h.m.n_params
		base_col = +~mod(j,2);
		base_row = row-2*floor((j-1)/2);
		h.sli_lab(j) = uicontrol('Style','Text','Units','normalized','String','INIT','Position',[0.05+base_col*0.5 row_offset(base_row-1)  0.4 rowheight],'Parent',h.fig,'HitTest','off');
		h.sli(j) = uicontrol('Style','slider','Units','normalized','Position',[0.05+base_col*0.5 row_offset(base_row)  0.4 rowheight],'Parent',h.fig,'Max',h.m.limits(2,j),'Min',h.m.limits(1,j),'Value',h.pars_orig(j));
		%set(h.sli(j),'Callback',@(a,b,c) redraw_callback(h.fig))
		h.listener(j) = addlistener(h.sli(j),'Value','PostSet',@(s,e)  redraw_callback(h.fig) );

	end
	guidata(h.fig,h)
	redraw_callback(h.fig)


	%hSlider = uicontrol('Style','slider','Callback',@(s,e) disp('hello'));



function redraw_callback(fig)
	h = guidata(fig);
	for j = 1:h.m.n_params
		set(h.sli_lab(j),'String',sprintf('%s = %.4f',h.m.param_names{j},get(h.sli(j),'Value')));
		pars(j) = get(h.sli(j),'Value');
	end
	draw(h.m,pars,h);


function draw(m,pars,h)

	try
		[f,P,stab] = m.spectrum(pars,m.target_f);
	catch
		f = NaN;
		P = NaN;
		stab = false;
	end
	a=get_pkf(f,P);
	%pks = [a.alpha_maxf,a.sigma_maxf,a.beta_maxf,a.gamma_maxf]

	beta_f = a.beta_maxf/a.alpha_maxf

	p = m.p_from_params(pars);
	alpha = p.alpha(1);
	beta = p.beta(1);
	t0 = p.t0;

	alpha_freq = 1/(t0 + 1/alpha + 1/beta); % Eq 21 HBM 2004
	%alpha_freq = 1/t0;
	sigma_freq = sqrt(beta*alpha)/2/pi;

	yl = get(h.ax,'YLim');
	set(h.s2,'XData',f,'YData',P);
	if stab
		set(h.s2,'LineStyle','--')
		set(h.tent_marker,'MarkerFaceColor','g','CData',[0 1 0])
	else
		set(h.s2,'LineStyle',':')
		set(h.tent_marker,'MarkerFaceColor','r','CData',[1 0 0])

	end

	set(h.t0_marker(1),'XData',[alpha_freq,alpha_freq],'YData',yl);
	set(h.t0_marker(2),'XData',[2*alpha_freq,2*alpha_freq],'YData',yl);
	set(h.t0_marker(3),'XData',[3*alpha_freq,3*alpha_freq],'YData',yl);
	set(h.spindle_marker,'XData',[sigma_freq,sigma_freq],'YData',yl);
	set(h.spindle_text,'Position',[sigma_freq,yl(1) 1]);
	set(h.t0_text(1),'Position',[alpha_freq,yl(1) 1]);
	set(h.t0_text(2),'Position',[2*alpha_freq,yl(1) 1]);
	set(h.t0_text(3),'Position',[3*alpha_freq,yl(1) 1]);

	[x,y,z,u] = tent.compute_mex(alpha,beta,t0,116);

	set(h.tent_surf,'XData',x,'YData',y,'ZData',z,'CData',u);
	set(h.tent_surf,'LineWidth',0.5,'FaceColor','none')
	edge_alphadata = x+y<0.90 & x > 0 & x < 1 & y > -1;
	set(h.tent_surf,'AlphaData',+edge_alphadata,'EdgeAlpha','flat','AlphaDataMapping','none')

	xyz = m.get_xyz(pars);
	set(h.tent_marker,'XData',xyz(1),'YData',xyz(2),'ZData',xyz(3));
	set(h.tent_drop,'XData',[xyz(1) xyz(1)],'YData',[xyz(2) xyz(2)],'ZData',[0 xyz(3)]);