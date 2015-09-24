function [cdata,states] = tent_cdata
	%ORIGINAL
	% cdata = [255 0 0; ...% EO bright red
	% 0 255 145; ...% EC light green
	% 255 180 0; ...% REM Orange
	% 0 159 100; ...% N1 dark green
	% 0 154 255; ...% N2 light blue
	% 12 23 178; ...% N3 dark blue
	% 102 0 51; ...% N2S purple
	% ];

	cdata = [137 44 145; ...% EO purple
	0 255 145; ...% EC light green
	255 180 0; ...% REM Orange
	255 0 0; ...% N1 Bright red
	0 154 255; ...% N2 light blue
	12 23 178; ...% N3 dark blue
	27 69 30; ...% N2S dark green
	];


	cdata = cdata/255;

	states = {'eo','ec','rem','n1','n2','n3','n2s'};