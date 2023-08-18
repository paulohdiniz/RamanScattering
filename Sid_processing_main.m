close all
clear variables 
cd('/Users/paulohd/Desktop/PauloDiniz')
addpath('/Users/paulohd/Desktop/PauloDiniz')

%% CODE FOR USING THE CLASS THAT MAKES THE SAME THING AS SID'S CODE

% We generate the class
SP=Sid_Processing();
SP.DCopt=1;

% Choose the experiment, and load all the data from the correct folders
SP=SP.choose_folders_load_data(3); %1,2,3

% Put the parameters
SP.ratio_window = 0.442;
SP.tukey_window_param = 0.5;
SP.deadtime=74;

% Removes the large curve before the sinusoidal
SP=SP.window_overlap_to_test(SP.tukey_window_param,SP.deadtime);

% Normalizes and centers the data
SP=SP.Tnorm_and_center_data(1,0,SP.deadtime);

% Chosen interpolation type
SP=SP.stitch_time_axis_T_with_interp('makima'); %makima, pchirp, linear

%
SP = SP.pick_fourier_window('blackman'); %blackman, tukeywin, hamming, hann, flattopwin, '...' for nothing

%FT
SP = SP.FT(SP.data_stitched.t_stitched, permute(SP.data_stitched.data_R,[3 1 2]).*repmat(SP.window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman spectrum is arbitrary units

SP = SP.make_raman_spectrum();

SP = SP.calculated_ssim_per_wn();

%Criteria ssym or frequency to find peaks
%SP = SP.points_to_plot_by_ssim();
SP = SP.points_to_plot_by_frequency();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(5 seconds)
plot_graphs(SP);

%(20 seconds)
%plot_graphs_with_filter(SP);

%(5 seconds)
%plot_graphs_with_mask(SP);

%(1 minute)
%x = get_best_ratio_window(SP, 0.2, 0.9);

%Save the txt with the results of filters (30 sec).
%get_best_filter(SP);

%(12 seconds)
%save_graphs_as_PDF(SP);