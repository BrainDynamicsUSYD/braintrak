function [cdata,state_names,short_mapping] = state_cdata

	cdata = [
	137  44  145; ...	% EO, purple
	0   255  145; ...	% EC, light  green
	255 180    0; ...	% REM, Orange
	255   0    0; ...	% N1, Bright red
	0   154  255; ...	% N2, light  blue
	12   23  178; ...	% N3, dark   blue
	27   69   30; ...	% N2S, dark   green
	0     0    0; ...   % Unresponsive, black
	];

	cdata = cdata/255;

	% These map directly to colours
	state_names = {'eo','ec','rem','n1','n2','n3','n2s','unresponsive'};

	% Map single character strings to colour indexes
	short_mapping{1}  = {'W','R','1','2','3','4'};
	short_mapping{2}  = [2 3 4 5 6 6]; % Indexes for cdata

	


	