function output = load_kaggle_data(set_type,set_number) 
    if nargin == 0
        set_type = 'preictal';
        set_number = 1;
    end

    window_length = 30;
    fft_length = 4;
    fname = sprintf('~/Desktop/psg_data/kaggle/patient_1_%s_%d.mat',set_type,set_number);
    fdata = load(fname);

    t = (0:(length(fdata.V)-1))/fdata.sampling_frequency;

    [output.t,output.f,output.s, output.nspec, output.n_reject] = bt.get_tfs(t,fdata.V,window_length,fft_length,false);
    
    output.colheaders = fdata.electrode;
    output.start_idx = 1;
    output.state_str = repmat({'W'},length(output.t),1);
    output.state_score = ones(size(output.state_str));
    output.Time = t;
    output.Cz = fdata.V;



	tic;
