classdef feather < handle
	% A feather is a container for a sequence of fit data
	% It contains an array of fit data and an array of plot data
	% corresponding to the outputs of bt.fit_spectrum()
	% It also contains a copy of the model used for the fit
	
	properties
		fit_data
		plot_data
		time
		path_prefix = ''
		model

		latest
		prealloc_size
        version
		viewer_obj
	end

	methods(Static)
		f = import(fname);
		%b = loadobj(a);
	end

	methods

		function self = feather(model,fit_data,plot_data,time)
			self.version = 3; % Keep track of changes to the interface
            
            if nargin < 4 || isempty(time)
            	time = 1;
            end

            if nargin == 0
            	return
            end

            self.model = model;

			self.fit_data = fit_data;
			if ischar(plot_data)
				[self.path_prefix,self.plot_data{1}] = fileparts(plot_data);
			else
				self.plot_data{1} = plot_data;
			end
			self.latest = 1;
			self.prealloc_size = 1;

			self.time = time;
		end

		function new = copy(self)
			error('This function is not completely functional')
			% The problem is that if there is a non-copyable property, it is
			% passed by reference. Specifically, this function will copy the feather
			% but the new feather's model will be a reference to the old one
			
			% Instantiate new object of the same class.
			new = feval(class(self));
			% Copy all non-hidden properties.
			p = properties(self);
			for i = 1:length(p)
				new.(p{i}) = self.(p{i});
			end
		end

		function append(self,f2)
			% Concatenate by adding another feather instance to the end
			assert(self.model == f2.model);
			for j = 1:f2.latest
				self.insert(f2.fit_data(j),f2.plot_data{j});
			end
		end

		function self = saveobj(self)
			% Automatically compress the object before saving
			self.compress();
		end

		function preallocate(self,n)
			self.prealloc_size = n; % Double the array size
			self.fit_data(n) = self.fit_data(1);
			self.plot_data{n} = self.plot_data{1};
		end

		function self = compress(self)
			self.prealloc_size = self.latest; % Double the array size
			self.fit_data = self.fit_data(1:self.latest);
			self.plot_data = self.plot_data(1:self.latest);
			self.time = self.time(1:self.latest);
		end

		function insert(self,fit_data,plot_data,time,idx)
			% Appends by default, overriding with an index lets you put things anywhere
			if self.latest == self.prealloc_size
				self.preallocate(self.latest*2);
			end

			if nargin < 4 || isempty(time)
				time = self.time(self.latest)+1; % Assume it is 1s after the last time
			end

			if nargin < 5 || isempty(idx) % Append the entry
				self.latest = self.latest + 1;
				idx = self.latest;
			elseif idx > self.latest + 1
				error('Cannot insert something that would leave a gap');
			end

			if ischar(plot_data)
				[test_prefix,self.plot_data{idx}] = fileparts(plot_data);
				assert(strcmp(test_prefix,self.path_prefix)==1);
			else
				self.plot_data{idx} = plot_data;
			end

			self.fit_data(idx) = fit_data;
			self.time(idx) = time;
		end

		function fill(self)
			% This function preloads all of the file 'pointers'
			% When f.plot_data = {'file_1','file_2'...}
			% f.fill will load the actual contents
			plot_data_new = {};
			for j = 1:length(self.plot_data)
				if ischar(self.plot_data{j})
					fdata = load(fullfile(self.path_prefix,self.plot_data{j}));
					plot_data_new{j} = fdata.plot_data;
				else
					plot_data_new{j} = self.plot_data{j};
				end
			end
			self.plot_data = plot_data_new;
		end


		function [fit_data,plot_data] = retrieve(self,idx)
			if idx > self.latest
				error('Index out of bounds');
			end
			fit_data = self.fit_data(idx);
            
            
			if ischar(self.plot_data{idx})
				fdata = load(fullfile(self.path_prefix,self.plot_data{idx}));
				plot_data = fdata.plot_data;
			else
				plot_data = self.plot_data{idx};
			end
		end

		function plot(self,idx)
			if (nargin < 2 || isempty(idx)) && self.latest > 1
				if isempty(self.viewer_obj) || self.viewer_obj.plot_chisq == 0 || ~ishandle(self.viewer_obj.parent_figure) 
					self.viewer_obj = bt.viewer(self);
				end
				return
			elseif (nargin < 2 || isempty(idx)) && self.latest == 1
				idx = 1;
			end

			[fd,pd] = self.retrieve(idx);
			f2 = bt.feather(self.model,fd,pd);
			if isempty(self.viewer_obj) || self.viewer_obj.plot_chisq == 1 || ~ishandle(self.viewer_obj.parent_figure) 
				self.viewer_obj = bt.viewer(f2);
			else
				self.viewer_obj.update_struct(fd,pd);
			end

		end

		function x = get_param(self,index)
			self.compress();
			x = zeros(self.latest,1);
			for j = 1:self.latest
				x(j) = self.fit_data(j).fitted_params(index);
			end
		end

	end
end

