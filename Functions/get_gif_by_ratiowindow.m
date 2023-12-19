function get_gif_by_ratiowindow(SPtemp)
    
    step = 0.02;
    ratios = 0.05:step:1;
    for i = 1:numel(ratios)
        SP(i) = SPtemp.copy(); %copying the data internally is much faster than loading all the files every time
        SP(i).ratio_window = ratios(i);
        SP(i) = SP(i).window_overlap_to_test(SP(i).tukey_window_param,SP(i).percent_FWHM);
        SP(i) = SP(i).Tnorm_and_center_data(1,0);
        SP(i) = SP(i).stitch_time_axis_T_with_interp(SP(i).interp_method);
        SP(i) = SP(i).pick_fourier_window(SP(i).window2_name); 
        SP(i) = SP(i).FT(SP(i).data_stitched.t_stitched, permute(SP(i).data_stitched.data_R,[3 1 2]).*repmat(SP(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
        SP(i) = SP(i).make_raman_spectrum();
        SP(i) = SP(i).points_to_plot_by_frequency();
        SP(i) = SP(i).get_signal_by_time_from_ifft();
        create_gif_with_ifft(SP(i), 'ratioWindow');    
        create_gif(SP(i), SP(i).window2_name);
        
    end

end


