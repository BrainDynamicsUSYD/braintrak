function animate_fit
	load steve
	load pp_allstates

	[time_bar,tent_point,tent_trail,spec_axis,tent_axis,exp_plot,fit_plot] = init_plot(t,f,s);

	xyz_trail = zeros(length(t),3);

	constraint_data = [];

	for j = 1:length(t)

		[f1,P1,idx] = db_fit(f,s(:,j),db_data,[]);
		keyboard
		title(spec_axis,sprintf('Step %d/%d',j,length(t)));
		xyz_trail(j,:) = db_data.xyz(idx,:);
		set(tent_point,'XData',xyz_trail(j,1),'YData',xyz_trail(j,2),'ZData',xyz_trail(j,3));
		set(tent_trail,'XData',xyz_trail(1:j,1),'YData',xyz_trail(1:j,2),'ZData',xyz_trail(1:j,3));

		set(exp_plot,'XData',log(f),'YData',log(s(:,j)));
		set(fit_plot,'XData',log(f1),'YData',log(P1));
		set(time_bar,'XData',[j j]);


		drawnow
		%fr = getframe(gcf);
		%writeVideo(vid,fr);
		constraint_data = update_constraints(db_data,idx);

	end

	%close(vid);
end

function constraint_data = update_constraints(db_data,fitted_points)
	% db_data is the points
	% fitted_points is an array with the index of all the points fitted so far

	% Based on db_data, use the current point to come up with weights 
	prev_xyz = mean(db_data.xyz(fitted_points,:));
	prev_nus = mean(db_data.nus(fitted_points,:));
	
	xyz_distance = bsxfun(@minus,db_data.xyz,prev_xyz(end,:)); 
	nus_distance = bsxfun(@minus,db_data.nus,prev_nus(end,:)); 

	xyz_fwhm = [0.15 0.4 0.2];
	nus_fwhm = 5e-3 * ones(1,8); % in parameter units
	to_std =  (2*sqrt(2*log(2)));

	constraint_data.xyz =  exp(-sum(bsxfun(@rdivide,xyz_distance.^2,2*(xyz_fwhm/to_std).^2),2));
	constraint_data.nus =  exp(-sum(bsxfun(@rdivide,nus_distance.^2,2*(nus_fwhm/to_std).^2),2));
end

function [time_bar,tent_point,tent_trail,spec_axis,tent_axis,exp_plot,fit_plot] = init_plot(t,f,s)
	fig = figure;
	%set(gcf,'Position',[821          56         421        1033])

	subplot(2,1,1)
	spec_axis = gca;
	%state_txt = text(0.2,-1.2,'TEST','FontSize',28,'FontWeight','bold');

	subplot(2,1,2);
	tent_axis = gca;
	analysis_to_tent(0,[],0.1,tent_axis);
	tent_point = scatter3(0,0,0,30,'ro','FaceColor','r');
	tent_trail = plot3(0,0,0,'r-');

	% subplot(3,1,3)
	% [tv,fv] = meshgrid(t,f);
	% [~,c] = contourf(tv,fv,log(s),50,'LineStyle','none');
	% %set(gca,'YScale','log')
	% colorbar
	% ylabel('Frequency (Hz)')
	% set(gcf,'Renderer','painters');
	% set(gca,'YScale','log');
	% set(gca,'YTick',[2:2:10 15:5:25 35 45]);
	% xl = get(gca,'XLim');
	% %set(gcf,'Position',[821          56         421        1033])
	% overlay = getframe;
	% delete(c);
	% imagesc([xl(1) xl(2)],[0 1],overlay.cdata);
	% set(gcf,'Renderer','opengl')
	% %set(gcf,'Position',[821          56         421        1033])
	% set(gca,'YTick',[])
	% hold on
	% time_bar = plot(gca,[0 0],[0 1],'r','LineWidth',3);
	time_bar = [];

	hold(spec_axis,'on');
	exp_plot = plot(spec_axis,log([1 2]),log([1 2]),'b');
	fit_plot = plot(spec_axis,log([1 2]),log([1 2]),'g');
	set(spec_axis,'XLim',log([1 45]),'XTick',log([1:8 10:2:16 20:5:45]),'XTickLabel',[1:8 10:2:16 20:5:45])
	set(spec_axis,'YTick',log(linspace(exp(-2),exp(4.5),20)),'YTickLabel','','YLim',[-2 4.5])
end