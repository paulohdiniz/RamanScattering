function [best_window, best_ratio_window, val_max] = get_best_window(SPtemp)

        val_temp = 0; 
        val_max = 0;
        best_indice = 1;
        windows = {'barthannwin', 'bartlett', 'blackman', 'blackmanharris', ...
                    'bohmanwin', 'chebwin', 'flattopwin', 'gausswin', 'hamming', ...
                    'hann', 'kaiser', 'nuttallwin', 'parzenwin', 'rectwin', ...
                    'taylorwin', 'tukeywin', 'triang', 'ones'};
        for i = 1:length(windows)
            SP(i) = SPtemp.copy(); %copying the data internally is much faster than loading all the files every time
            SP(i).window2_name = windows{i};
            SP(i).ratio_window = get_best_ratio_window_by_frequency(SP(i));
            SP(i) = SP(i).window_overlap_to_test(SP(i).tukey_window_param,SP(i).deadtime);
            SP(i) = SP(i).Tnorm_and_center_data(1,0,SP(i).deadtime);
            SP(i) = SP(i).stitch_time_axis_T_with_interp(SP(i).interp_method);
            SP(i) = SP(i).pick_fourier_window(SP(i).window2_name); 
            SP(i) = SP(i).FT(SP(i).data_stitched.t_stitched, permute(SP(i).data_stitched.data_R,[3 1 2]).*repmat(SP(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
            SP(i) = SP(i).make_raman_spectrum();
            SP(i) = SP(i).points_to_plot_by_frequency();
            SP(i).scoreCriteria = getScoreCriteria(SP(i));
            
            val_temp = SP(i).scoreCriteria(2);
            
            if(val_max < val_temp)
                val_max = val_temp;
                best_indice = i;
            end
            
            best_window = SP(best_indice).window2_name;
            best_ratio_window = SP(best_indice).ratio_window;
        end

end

