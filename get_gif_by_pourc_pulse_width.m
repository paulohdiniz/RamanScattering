function get_gif_by_pourc_pulse_width(SPtemp)
    step = 2;
    ratios = 50:step:150;
    for i = 1:numel(ratios)
        SP(i) = SPtemp.copy(); %copying the data internally is much faster than loading all the files every time
        SP(i).pourc_pulse_width = ratios(i);
        SP(i) = SP(i).window_overlap_to_test(SP(i).tukey_window_param,SP(i).pourc_pulse_width);
        SP(i) = SP(i).Tnorm_and_center_data(1,0);
        SP(i) = SP(i).stitch_time_axis_T_with_interp(SP(i).interp_method);
        SP(i) = SP(i).pick_fourier_window(SP(i).window2_name); 
        SP(i) = SP(i).FT(SP(i).data_stitched.t_stitched, permute(SP(i).data_stitched.data_R,[3 1 2]).*repmat(SP(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
        SP(i) = SP(i).make_raman_spectrum();
        SP(i) = SP(i).points_to_plot_by_frequency();
        SP(i).signalIFFT = get_time_spec_from_peaks_by_ifft(SP(i));
        create_gif_with_ifft(SP(i), 'pourc_pulse');    
        %create_gif(SP(i), 'Deadtime');
    end
   
end


