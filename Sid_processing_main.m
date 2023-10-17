close all
clear variables 
%cd('/Users/paulohd/Desktop/PauloDiniz')
cd('C:\Users\phdin\Desktop\Data_Paulo\210507') % '210306' (xp 29, 30), 2ele
                                               % '210507' (xp 41, 15, 16, 17) 2ele
                                               % '210617' (xp 4, 19), 1 ele
%addpath('/Users/paulohd/Desktop/PauloDiniz')  %210619 xp 10 pour analyser
addpath('C:\Users\phdin\Desktop\Data_Paulo\210507')

%% CODE FOR USING THE CLASS THAT MAKES THE SAME THING AS SID'S CODE

% We generate the class
SP=Sid_Processing();
SP.DCopt=1;

% Choose the experiment, and load all the data from the correct folders
SP=SP.choose_folders_load_data(41); 

% Put the parameters
SP.tukey_window_param = 1;
SP.deadtime=70; %110
SP.window2_name = 'rectwin'; %barthannwin, bartlett, blackman, blackmanharris, bohmanwin, 
                        % chebwin, flattopwin, gausswin, hamming, hann,
                        % kaiser, nuttallwin, parzenwin, rectwin,
                        % taylorwin, tukeywin,tukeywinINV, triang, ones
SP.interp_method = 'makima'; %makima, pchirp, spline

%[SP.window2_name, SP.ratio_window] = get_best_window(SP);
SP.ratio_window = 1;%get_best_ratio_window_by_frequency(SP); % after tukey and deadtime
%SP.ratio_window = get_best_ratio_window_by_ssim(SP); % after tukey and deadtime

% Removes the large curve before the sinusoidal
SP=SP.window_overlap_to_test(SP.tukey_window_param,SP.deadtime);
%SP=SP.window_overlap(SP.tukey_window_param,SP.deadtime); %change de 

% Normalizes and centers the data
SP=SP.Tnorm_and_center_data(1,0,SP.deadtime);

% Chosen interpolation type
SP=SP.stitch_time_axis_T_with_interp(SP.interp_method);

%
SP = SP.pick_fourier_window(SP.window2_name); 

%FT
SP = SP.FT(SP.data_stitched.t_stitched, permute(SP.data_stitched.data_R,[3 1 2]).*repmat(SP.window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman spectrum is arbitrary units

SP = SP.make_raman_spectrum();

SP = SP.calculated_images_scores_per_wn();

%Criteria ssym or frequency to find peaks

SP = SP.points_to_plot_by_frequency();
%SP = SP.points_to_plot_by_ssim();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get_gif_by_deadtime(SP);
get_gif_by_interp_method(SP);
%get_gif_by_ratiotukey(SP);


%get_gif_by_window(SP);
%get_model(SP);
%get_time_spec_from_peaks_by_ifft(SP);


%PLOT 
%get_time_spec_from_peaks_by_ifft(SP);
%plot_graphs(SP); %(5 seconds)
%plot_graphs_with_transmission(SP);
%plot_graphs_with_filter(SP); %(20 seconds)
%plot_graphs_with_mask(SP); %(5 seconds)
%plot_best_ssim_by_ratio_window(SP);
%plot_spectrogram(SP);
%save_graphs_as_PDF(SP); %(12 seconds)
%save_raman_spectrum_as_PDF(SP);
%save_windows_by_ratio_as_gif(SP);



%RATIO WINDOW
%x = get_best_ratio_window_by_ssim(SP); %(1 min)
%y = get_best_ratio_window_by_frequency(SP); %(10 sec)

%Save the txt with the results of filters (30 sec).
%get_best_filter(SP);


