function best_ratio = get_best_ratio_window_by_frequency(SPtemp)

%analyse ratiowindow
    step = 0.05;
    ratios = 0.05:step:1;
    for i = 1:numel(ratios)
        SP(i) = SPtemp.copy(); %copying the data internally is much faster than loading all the files every time
        SP(i).ratio_window = ratios(i);
        SP(i) = SP(i).window_overlap_to_test(SP(i).tukey_window_param,SP(i).deadtime);
        SP(i) = SP(i).Tnorm_and_center_data(1,0,SP(i).deadtime);
        SP(i) = SP(i).stitch_time_axis_T_with_interp(SP(i).interp_method);
        SP(i) = SP(i).pick_fourier_window(SP(i).window2_name); 
        SP(i) = SP(i).FT(SP(i).data_stitched.t_stitched, permute(SP(i).data_stitched.data_R,[3 1 2]).*repmat(SP(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
        SP(i) = SP(i).make_raman_spectrum();
        SP(i) = SP(i).points_to_plot_by_frequency();
        
    end
    
    best_indice = 1; %wanted interaction 
    val_temp = 0.0; 
    val_max = 0.0;
    for i = 1:numel(ratios)

       val_temp = SP(i).scoreCriteria(2) * ~isempty(SP(i).scoreCriteria) + ...
            isempty(SP(i).scoreCriteria) * 0.0;
        
       if(val_max < val_temp)
            val_max = val_temp;
            best_indice = i;
       end
        
    end
    best_ratio_temp = SP(best_indice).ratio_window;
    
    % Calcul de ratio_init et ratio_fin
    ratio_init = (best_ratio_temp ~= 0.05) * (best_ratio_temp - 0.05) + (best_ratio_temp == 0.05) * 0.5;
    ratio_fin = (best_ratio_temp ~= 1) * (best_ratio_temp + 0.05) + (best_ratio_temp == 1) * 1;

    step = 0.01;

    ratios = ratio_init:step:ratio_fin;
    for i = 1:numel(ratios)
        SP(i) = SPtemp.copy(); %copying the data internally is much faster than loading all the files every time
        SP(i).ratio_window = ratios(i);
        SP(i) = SP(i).window_overlap_to_test(SP(i).tukey_window_param,SP(i).deadtime);
        SP(i) = SP(i).Tnorm_and_center_data(1,0,SP(i).deadtime);
        SP(i) = SP(i).stitch_time_axis_T_with_interp(SP(i).interp_method); %makima, pchirp, linear
        SP(i) = SP(i).pick_fourier_window(SP(i).window2_name); %SP(i).window2_name , blackman, tukeywin, hamming, hann, flattopwin, '...' for nothing
        SP(i) = SP(i).FT(SP(i).data_stitched.t_stitched, permute(SP(i).data_stitched.data_R,[3 1 2]).*repmat(SP(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
        SP(i) = SP(i).make_raman_spectrum();
        SP(i) = SP(i).points_to_plot_by_frequency();
        
    end
    
    best_indice = 1; %wanted interaction 
    val_temp = 0.0; 
    val_max = 0.0;
    for i = 1:numel(ratios)

       val_temp = SP(i).scoreCriteria(2) * ~isempty(SP(i).scoreCriteria) + ...
            isempty(SP(i).scoreCriteria) * 0.0;
        
       if(val_max < val_temp)
            val_max = val_temp;
            best_indice = i;
       end
        
    end
    best_ratio_temp = SP(best_indice).ratio_window;

    fprintf('The best ratio window(%s) by frequency is %s.\n',string(SP(best_indice).window2_name),string(best_ratio_temp))
    best_ratio = best_ratio_temp;
    
end

