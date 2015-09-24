function [x,y,theta,l] = electrode_positions(electrodes,Lx,Ly)
	% Retrieve the electrode positions given a cell array of electrode_positions

	% Note that l=0.5 needs to correspond to 0.1995 units of distance
	dscale = 0.25; % Original scaling
	%dscale = 0.3990; % Circle of radius 0.1995m has area equal to half the 50x50cm grid

	if nargin == 0 && nargout == 1
		x = dscale;
		return
	end

	if nargin < 3 || isempty(Ly)
		Ly = 0.5;
	end

	if nargin < 2 || isempty(Lx)
		Lx = 0.5;
	end

	if isstr(electrodes)
		electrodes = {electrodes};
	end
	
	x = zeros(length(electrodes),1);
	y = x;
	theta = x;
	l = x;

	for j = 1:length(electrodes)
		[x(j),y(j),theta(j),l(j)] = get_eeglab_position(upper(electrodes{j}),Lx,Ly,dscale);
	end

function [x,y,theta,l] = get_eeglab_position(electrode,Lx,Ly,dscale)
	% Drawn from EEGLAB
	% Head circumference = 2 for these data
	% Thus an arc length (second column) of 0.5 corresponds to going 
	% to the side or front of the head. Angle (first column) of 0
	% corresponds to the forward direction
	electrode = upper(electrode);

	switch electrode
		case 'T3'
			electrode = 'T7';
		case 'T4'
			electrode = 'T8';
		case 'T5'
			electrode = 'P7';
		case 'T6'
			electrode = 'P8';

	end

	data = {
		0       , 0.50669 , 'FPZ'
		23      , 0.71    , 'EOG1'
		-39.947 , 0.34459 , 'F3'
		0       , 0.25338 , 'FZ'
		39.897  , 0.3445  , 'F4'
		-43     , 0.65    , 'EOG2'
		-69.332 , 0.40823 , 'FC5'
		-44.925 , 0.18118 , 'FC1'
		44.925  , 0.18118 , 'FC2'
		69.332  , 0.40823 , 'FC6'
		-90     , 0.53318 , 'T7'
		-90     , 0.26669 , 'C3'
		90      , 0.26667 , 'C4'
		0       , 0       , 'CZ'
		90      , 0.53318 , 'T8'
		-110.67 , 0.40823 , 'CP5'
		-135.07 , 0.18118 , 'CP1'
		135.07  , 0.18118 , 'CP2'
		110.67  , 0.40823 , 'CP6'
		-126.09 , 0.52808 , 'P7'
		-140.05 , 0.34459 , 'P3'
		180     , 0.25338 , 'PZ'
		140.1   , 0.3445  , 'P4'
		126.13  , 0.52807 , 'P8'
		-144.11 , 0.52233 , 'PO7'
		-157.54 , 0.42113 , 'PO3'
		180     , 0.37994 , 'POZ'
		157.54  , 0.42113 , 'PO4'
		144.14  , 0.52231 , 'PO8'
		-162.07 , 0.51499 , 'O1'
		180     , 0.50669 , 'OZ'
		162.07  , 0.51499 , 'O2'
		0       , 0.4     , 'TEST_FRONT'
		0       , 0       , 'TEST_MID'
		180     , 0.4     , 'TEST_REAR'
		};

		% From ftp://sccn.ucsd.edu/pub/locfiles/eeglab/chan61fms11.loc
		data = { 
			-18    , 0.415 , 'FP1'
			0      , 0.4   , 'FPZ'
			18     , 0.415 , 'FP2'
			-35    , 0.405 , 'AF7'
			-22    , 0.32  , 'AF3'
			0      , 0.3   , 'AFZ'
			22     , 0.32  , 'AF4'
			35     , 0.405 , 'AF8'
			-53.5  , 0.412 , 'F7'
			-47    , 0.338 , 'F5'
			-38    , 0.27  , 'F3'
			-21    , 0.219 , 'F1'
			0      , 0.2   , 'FZ'
			21     , 0.219 , 'F2'
			38     , 0.27  , 'F4'
			47     , 0.338 , 'F6'
			53.5   , 0.412 , 'F8'
			-70    , 0.4   , 'FT7'
			-67    , 0.3   , 'FC5'
			-59    , 0.215 , 'FC3'
			-41    , 0.138 , 'FC1'
			0      , 0.1   , 'FCZ'
			41     , 0.138 , 'FC2'
			59     , 0.215 , 'FC4'
			67     , 0.3   , 'FC6'
			70     , 0.4   , 'FT8'
			-90    , 0.4   , 'T7'
			-90    , 0.3   , 'C5'
			-90    , 0.2   , 'C3'
			-90    , 0.1   , 'C1'
			-90    , 0     , 'CZ'
			90     , 0.1   , 'C2'
			90     , 0.2   , 'C4'
			90     , 0.3   , 'C6'
			90     , 0.4   , 'T8'
			-110   , 0.4   , 'TP7'
			-113   , 0.3   , 'CP5'
			-121   , 0.215 , 'CP3'
			180    , 0.1   , 'CPZ'
			-139   , 0.138 , 'CP1'
			139    , 0.138 , 'CP2'
			121    , 0.215 , 'CP4'
			113    , 0.3   , 'CP6'
			110    , 0.4   , 'TP8'
			-126.5 , 0.412 , 'P7'
			-133   , 0.338 , 'P5'
			-142   , 0.27  , 'P3'
			-159   , 0.219 , 'P1'
			180    , 0.2   , 'PZ'
			159    , 0.219 , 'P2'
			142    , 0.27  , 'P4'
			133    , 0.338 , 'P6'
			126.5  , 0.412 , 'P8'
			-145   , 0.405 , 'PO7'
			-158   , 0.32  , 'PO3'
			180    , 0.3   , 'POZ'
			158    , 0.32  , 'PO4'
			145    , 0.405 , 'PO8'
			-162   , 0.415 , 'O1'
			180    , 0.4   , 'OZ'
			162    , 0.415 , 'O2'
			0       , 0.5     , 'TEST_FRONT'
			0       , 0       , 'TEST_MID'
			180     , 0.5     , 'TEST_REAR'
		};

	% Convert the EEGLAB position to X and Y coordinates
	idx = find(strcmp(electrode,data(:,3)));
	theta = data{idx,1};
	l = data{idx,2};


	x = Lx/2 + l*dscale*cosd(-theta+90);
	y = Ly/2 - l*dscale*sind(-theta+90); % Because in the model framework, y=0 is the front