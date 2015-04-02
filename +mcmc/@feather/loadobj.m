function loaded = loadobj(loaded)
    % This function is designed to convert older instances
    % into newer instances
    return

   	% Error is: "You cannot set the read-only property 'param_symbols' of full."
    if isempty(loaded.model.param_symbols)
    	disp('Upgrading model parameter symbols and units using same class')
    	m2 = feval(class(loaded.model));
    	loaded.model.param_symbols = m2.param_symbols;
    	loaded.model.param_units = m2.param_units;
    end

end

