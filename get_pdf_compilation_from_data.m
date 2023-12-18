close all
clear variables 
cd('/Users/paulohd/Documents/Data_SID/') % '210306' (xp 29, 30), 2ele
                                         % '210507' (xp 41, 15, 16, 17) 2ele
                                         % '210617' (xp 4, 19), 1 ele
addpath('/Users/paulohd/Documents/Data_SID/')
addpath('/Users/paulohd/Desktop/PauloDiniz')
%addpath(genpath('/Users/paulohd/Documents/Data_SID/'))
%%%%%%%%%
Data_path = '/Users/paulohd/Documents/Data_SID/';
folder_content_dates=dir(Data_path); % We look a the folder content
folder_content_dates = folder_content_dates(~startsWith({folder_content_dates.name}, ".")); %excludes all starting with "."
folder_content_dates = folder_content_dates([folder_content_dates.isdir]); %only folders, not files.


%75 elements
for i = 1:numel(folder_content_dates)
    cd(append(folder_content_dates(i).folder, '/', folder_content_dates(i).name))
    folder_content = dir(cd);
    folder_content = folder_content(~startsWith({folder_content.name}, ".")); %excludes all starting with "."
    folder_content = folder_content([folder_content.isdir]); %only folders, not files.

    names_for_dir=cellfun(@(x) x(1:numel(x)),{folder_content(:).name},'UniformOutput',false);  
    % names=cellfun(@(x) x(9:145),{folder_content(:).name},'UniformOutput',false);  

    % We generate the class
    SP=Sid_Processing();
    SP.name_pdf = folder_content_dates(i).name;

    % Choose the experiment, and load all the data from the correct folders
    for k=1:numel(names_for_dir)/5
        SP=SP.choose_folders_load_data(k); %1,2,3,..
        %SP.ratio_window = get_best_ratio_window_by_frequency(SP);
        SP.ratio_window = 0.5;
        SP.tukey_window_param = 1;
        SP.pourc_pulse_width=100;
        SP=SP.window_overlap_to_test(SP.tukey_window_param,SP.percent_FWHM);
        SP=SP.Tnorm_and_center_data(1,0);
        SP=SP.stitch_time_axis_T_with_interp('makima'); %makima, pchirp, linear
        SP = SP.pick_fourier_window('blackman'); %blackman, tukeywin, hamming, hann, flattopwin, '...' for nothing
        SP = SP.FT(SP.data_stitched.t_stitched, permute(SP.data_stitched.data_R,[3 1 2]).*repmat(SP.window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman spectrum is arbitrary units
        SP = SP.make_raman_spectrum();
        SP = SP.calculated_images_scores_per_wn();
        SP = SP.points_to_plot_by_frequency();
        save_graphs_as_PDF(SP);
    end
end


%desse jeito dentro de cada experimento alem do 0000.txt e outros arquivos,
%teremos um arquivo texto com nome especial com todos os parametros do
%titutlo.

%feito isso, podemos entrar em cada experimento rodar o Sid_processing_main
%e jogar a saida do matlab num pdf, mas antes dizendo os parametros que
%achamos no titulo.

