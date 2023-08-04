close all
clear variables 
cd('/Users/paulohd/Documents/Data_SID/')
addpath('/Users/paulohd/Documents/Data_SID/')
%%%%%%%%%
Data_path = '/Users/paulohd/Documents/Data_SID/';
folder_content_dates=dir(Data_path); % We look a the folder content
folder_content_dates = folder_content_dates(~startsWith({folder_content_dates.name}, ".")); %excludes all starting with "."
folder_content_dates = folder_content_dates([folder_content_dates.isdir]); %only folders, not files.


%laço pra entrar em todos, vou fazer com o 150 aqui
cd(append(folder_content_dates(150).folder, '/', folder_content_dates(150).name))
folder_content = dir(cd);
folder_content = folder_content(~startsWith({folder_content.name}, ".")); %excludes all starting with "."
folder_content = folder_content([folder_content.isdir]); %only folders, not files.

%salvei em foldercontent todas as pastas dentro
names_for_dir=cellfun(@(x) x(1:numel(x)),{folder_content(:).name},'UniformOutput',false);  
names=cellfun(@(x) x(9:145),{folder_content(:).name},'UniformOutput',false);  
%em names salvei o nome dos foldes

for i=1:numel(names)
end

psdelay = string(regexp(names_for_dir{1}, '_([0-9]+)psdelay', 'tokens')); % match prendre tout
function_generator = string(regexp(names_for_dir{1}, 'FG_([0-9]+)MHZ', 'match'));
lockin_parameters = string(regexp(names_for_dir{1}, 'APE([\w]+)ns', 'match'));
%objectives = regexp(names{1}, 'TODO', 'match'); REFAIRE

oldfolder = cd (append(folder_content(1).folder, '/' , names_for_dir{1}));
fileID = fopen('parameters_title.txt','w');
fprintf(fileID, 'Date: %s \n', datetime("now"));
fprintf(fileID, 'Intern: Paulo Henrique DINIZ FERNANDES\n');
fprintf(fileID, 'psdelay: %s \n', psdelay);
fprintf(fileID, 'Function_generator: %s \n', function_generator);
fprintf(fileID, 'Lock-in_parameters: %s \n', lockin_parameters);
fprintf(fileID, '-------------------------------------------\n');
fclose(fileID);
cd(oldfolder)


% We generate the class
SP=Sid_Processing();
SP.DCopt=1;
% Choose the experiment, and load all the data from the correct folders
for i=1:numel(names)/5
    i
    tic
    SP=SP.choose_folders_load_data(i); %1,2,3,..
    toc
    SP.ratio_window = get_best_ratio_window(SP, 0.1, 0.9);
    SP.tukey_window_param = 1;
    SP.deadtime=74;
    SP=SP.window_overlap_to_test(SP.tukey_window_param,SP.deadtime);
    SP=SP.Tnorm_and_center_data(1,0,SP.deadtime);
    SP=SP.stitch_time_axis_T_with_interp('makima'); %makima, pchirp, linear
    SP = SP.pick_fourier_window('blackman'); %blackman, tukeywin, hamming, hann, flattopwin, '...' for nothing
    SP = SP.FT(SP.data_stitched.t_stitched, permute(SP.data_stitched.data_R,[3 1 2]).*repmat(SP.window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman spectrum is arbitrary units
    SP = SP.make_raman_spectrum();
    SP = SP.points_to_plot_by_ssim();
    save_graphs_as_PDF(SP);
end


%desse jeito dentro de cada experimento alem do 0000.txt e outros arquivos,
%teremos um arquivo texto com nome especial com todos os parametros do
%titutlo.

%feito isso, podemos entrar em cada experimento rodar o Sid_processing_main
%e jogar a saida do matlab num pdf, mas antes dizendo os parametros que
%achamos no titulo.

