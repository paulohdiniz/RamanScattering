close all
clear variables 
cd('C:\Users\phdin\Desktop\Data_Paulo\210507\') % '210306' (xp 29, 30), 2ele
                                               % '210507' (xp 41, 15, 16, 17) 2ele
                                               % '210617' (xp 4, 19), 1 ele
%addpath('/Users/paulohd/Desktop/210619(1samp)')  %210619 xp 10 pour analyser
addpath('C:\Users\phdin\Desktop\PauloDiniz')

%% Script that manages the Raman experiment

% We generate the class
XP = Raman_Processing();

% Choose the experiment, and load all the data from the correct folders
XP=XP.choose_folders_load_data(15);

RP=RmFast(XP.data_raw, XP.N_x, XP.N_y, XP.N_t, XP.Clock_Freq);

% % Put the parameters
% RP.tukey_window_param = 1;
% RP.pourc_pulse_width=100;
% 
% % Second window
% RP.window2_name = 'blackman'; %barthannwin, bartlett, blackman, blackmanharris, bohmanwin,
%                         % chebwin, flattopwin, gausswin, hamming, hann,
%                         % kaiser, nuttallwin, parzenwin, rectwin, 
%                         % taylorwin, tukeywin,tukeywinINV, triang, ones
%
%
% %Ratio of second window
% RP.ratio_window = 1;
% 
% % Removes the large curve before the sinusoidal
% RP=RP.window_overlap(RP.tukey_window_param,RP.pourc_pulse_width);
% 
% % Normalizes and centers the data
% RP=RP.Tnorm_and_center_data();
% 
% % Do the interpolation
% RP=RP.stitch_time_axis_T_with_interp();
% 
% % Applies the second window to the signal
% RP = RP.pick_fourier_window(RP.window2_name);
% 
% % FT
% RP = RP.FT(RP.data_stitched.t_stitched, permute(RP.data_stitched.data_R,[3 1 2]).*repmat(RP.window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman spectrum is arbitrary units
% 
% RP = RP.make_raman_spectrum();

RP = RP.points_to_plot_by_frequency();

imagesc(squeeze(abs(RP.hyperspectralRamanImageComplex(RP.pixels_plot(1),:,:))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GIFS
%get_gif_by_pourc_pulse_width(RP);
%get_gif_by_interp_method(RP);
%get_gif_by_ratiotukey(RP);
%get_gif_by_ratiowindow(RP);
%get_gif_by_windows(RP);
%get_infos_by_windows(RP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT
%plot_graphs(RP);
%plot_graphs_with_transmission(RP);
%plot_graphs_with_roi(RP);
plot_graphs_with_ifft(RP);
%plot_best_ssim_by_ratio_window(RP);
%plot_spectrogram(RP);
%plot_best_hyperspectral_images(RP);
%plot_images_with_filters_by_psnr(RP);
%plot_similar_pixels_from_rois(RP);

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