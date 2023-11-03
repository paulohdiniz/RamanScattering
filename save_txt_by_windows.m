function save_txt_by_windows(SPtemp)
    
    fileID = fopen(append('infos_windows_exp',string(SPtemp.xp_number), '.txt'),'w');
    fprintf(fileID, 'Analysis date: %s \n', datetime("now"));
    fprintf(fileID, 'Intern: Paulo Henrique DINIZ FERNANDES\n');
    fprintf(fileID, 'Folder path : %s \n', string(pwd));
    fprintf(fileID, 'Experiment number: %d \n', SPtemp.xp_number);
    fprintf(fileID, 'Tukey Window parameter: %.4f \n', SPtemp.tukey_window_param);
    fprintf(fileID, 'Pourc pulse width: %.4f \n', SPtemp.pourc_pulse_width);
    fprintf(fileID, '\n');
    windows = {'barthannwin', 'bartlett', 'blackman', 'blackmanharris', ...
                'bohmanwin', 'chebwin', 'flattopwin', 'gausswin', 'hamming', ...
                'hann', 'kaiser', 'nuttallwin', 'parzenwin', 'rectwin', ...
                'taylorwin', 'tukeywin', 'triang', 'ones'};

    for i = 1:length(windows)
        SP(i) = SPtemp.copy(); %copying the data internally is much faster than loading all the files every time
        SP(i).window2_name = windows{i};
        SP(i) = SP(i).window_overlap_to_test(SP(i).tukey_window_param,SP(i).pourc_pulse_width);
        SP(i) = SP(i).Tnorm_and_center_data(1,0);
        SP(i) = SP(i).stitch_time_axis_T_with_interp(SP(i).interp_method);
        SP(i) = SP(i).pick_fourier_window(SP(i).window2_name); 
        SP(i) = SP(i).FT(SP(i).data_stitched.t_stitched, permute(SP(i).data_stitched.data_R,[3 1 2]).*repmat(SP(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
        SP(i) = SP(i).make_raman_spectrum();
        SP(i) = SP(i).points_to_plot_by_frequency();
        SP(i).signalIFFT = get_time_spec_from_peaks_by_ifft(SP(i));


        %TXT WITH RESULTS

        fprintf(fileID, 'Window type: %s \n', SP(i).window2_name);
        fprintf(fileID, 'Ratio Window: %.2f \n', SP(i).ratio_window);
        fprintf(fileID, 'Peaks in order: %d, %d, %d \n', SP(i).pixels_plot(1), SP(i).pixels_plot(2), SP(i).pixels_plot(3));
        fprintf(fileID, 'Peaks in order (cm-1): %.2f, %.2f, %.2f \n', SP(i).wn(SP(i).pixels_plot(1)), SP(i).wn(SP(i).pixels_plot(2)), SP(i).wn(SP(i).pixels_plot(3)));
        
        fprintf(fileID, 'Peak 1: %.3f, %.2e, %.2e, %.2f, %.2e\n', SP(i).IP.ssim_wn(SP(i).pixels_plot(1)),SP(i).peakProm(1), SP(i).peakAmpli(1), SP(i).peakWidth(1), SP(i).peakProm(1)/SP(i).peakWidth(1));
        fprintf(fileID, 'Peak 2: %.3f, %.2e, %.2e, %.2f, %.2e\n', SP(i).IP.ssim_wn(SP(i).pixels_plot(2)),SP(i).peakProm(2), SP(i).peakAmpli(2), SP(i).peakWidth(2), SP(i).peakProm(2)/SP(i).peakWidth(2));
        fprintf(fileID, 'Peak 3: %.3f, %.2e, %.2e, %.2f, %.2e\n', SP(i).IP.ssim_wn(SP(i).pixels_plot(3)),SP(i).peakProm(3), SP(i).peakAmpli(3), SP(i).peakWidth(3), SP(i).peakProm(3)/SP(i).peakWidth(3));
        fprintf(fileID, '\n');
        %create_gif_with_ifft(SP(i), 'windows');    
        %create_gif(SP(i), 'windows');
    end
        fclose(fileID);  
end
