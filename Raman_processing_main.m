close all
clear variables 
cd('C:\Users\phdin\Desktop\Data_Paulo\210507') % '210306' (xp 29, 30), 2ele
                                               % '210507' (xp 41, 15, 16, 17) 2ele
                                               % '210617' (xp 4, 19), 1 ele
%addpath('/Users/paulohd/Desktop/210619(1samp)')  %210619 xp 10 pour analyser
addpath('C:\Users\phdin\Desktop\Data_Paulo\210507')
xp_number = 41;

%% Script that manages the Raman experiment

% We generate the class
RP=Raman_Processing();

% Choose the experiment, and load all the data from the correct folders
RP=RP.choose_folders_load_data(xp_number); 

% Put the parameters
RP.tukey_window_param = 1;
RP.percent_FWHM=100;

% Interpolation method
RP.interp_method = 'makima'; % nearest, linear, makima, pchip, spline

% Second window
RP.window2_name = 'blackman'; %barthannwin, bartlett, blackman, blackmanharris, bohmanwin, 
                        % chebwin, flattopwin, gausswin, hamming, hann,
                        % kaiser, nuttallwin, parzenwin, rectwin,
                        % taylorwin, tukeywin,tukeywinINV, triang, ones

%Ratio of second window
%[RP.window2_name, RP.ratio_window] = get_best_window(RP);
RP.ratio_window = 1;%get_best_ratio_window_by_frequency(RP); % after tukey and deadtime
%RP.ratio_window = get_best_ratio_window_by_ssim(RP); % after tukey and deadtime

% Removes the large curve before the sinusoidal
RP=RP.window_overlap_to_test(RP.tukey_window_param,RP.percent_FWHM);
%RP=RP.window_overlap(RP.tukey_window_param,RP.percent_FWHM); %change de 

% Normalizes and centers the data
RP=RP.Tnorm_and_center_data(1,0);

% Do the interpolation
RP=RP.stitch_time_axis_T_with_interp(RP.interp_method);

% Applies the second window to the signal
RP = RP.pick_fourier_window(RP.window2_name); 

% Creating hyperspectral datacube with Fourier Transform
RP = RP.FT(RP.data_stitched.t_stitched, permute(RP.data_stitched.data_R,[3 1 2]).*repmat(RP.window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman spectrum is arbitrary units

% Raman spectrum from hyperspectral datacube
RP = RP.make_raman_spectrum();

% Calculates image quality indices
RP = RP.calculated_images_scores_per_wn();

% Find peaks by Prominence
RP = RP.points_to_plot_by_frequency();

% Creates signal by time from ifft of peaks
RP = RP.get_signal_by_time_from_ifft();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT 
% plot_graphs(RP);
% plot_graphs_with_transmission(RP);
% plot_graphs_with_roi(RP);
% plot_similar_pixels_from_rois(RP);
% plot_graphs_with_ifft(RP);
% plot_best_ssim_by_ratio_window(RP);
% plot_spectrogram(RP);
% plot_best_hyperspectral_images(RP);
% plot_images_with_filters_by_psnr(RP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GIFS
get_gif_by_FWHM(RP);
get_gif_by_interp_method(RP);
get_gif_by_ratiotukey(RP);
get_gif_by_ratiowindow(RP);
get_gif_by_windows(RP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE
%save_graphs_as_PDF(RP); 
%save_raman_spectrum_as_PDF(RP);
%save_txt_by_windows(RP);
%save_txt_best_filter_by_ssim(RP);
%save_txt_best_filter_by_brisque(RP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TO DO:
%get_model(RP);




