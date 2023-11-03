function best_ratio = get_best_ratio_window_by_ssim(SPtemp)

%analyse ratiowindow
    step = 0.05;
    ratios = 0.05:step:1;
    for i = 1:numel(ratios)
        SP(i) = SPtemp.copy(); %copying the data internally is much faster than loading all the files every time
        SP(i).ratio_window = ratios(i);
        SP(i) = SP(i).window_overlap_to_test(SP(i).tukey_window_param,SP(i).pourc_pulse_width);
        SP(i) = SP(i).Tnorm_and_center_data(1,0);
        SP(i) = SP(i).stitch_time_axis_T_with_interp(SP(i).interp_method); %makima, pchirp, linear
        SP(i) = SP(i).pick_fourier_window(SP(i).window2_name); %SP(i).window2_name , blackman, tukeywin, hamming, hann, flattopwin, '...' for nothing
        SP(i) = SP(i).FT(SP(i).data_stitched.t_stitched, permute(SP(i).data_stitched.data_R,[3 1 2]).*repmat(SP(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
        SP(i) = SP(i).make_raman_spectrum();
    end
    
    best_indice = 1; %wanted interaction 
    val_temp = 0; 
    val_max = 0;
    for i = 1:numel(ratios)
        SP(i).IP = Image_Processing();

        temp2=squeeze(mean((SP(i).data_processed(1).data_T),[1]));
        temp2=(temp2-min(temp2(:)))./max(temp2(:));
        SP(i).IP.mat_ref = temp2;

        for j = 1:size(SP(i).hyperspectralRamanImageComplex,1)/2
            temp=squeeze(abs(SP(i).hyperspectralRamanImageComplex(j,:,:)));
            temp=(temp-min(temp(:)))./max(temp(:));
            SP(i).IP.mat_img_wn{j} = temp;
            SP(i).IP.ssim_wn(j) = ssim(temp,temp2);
        end
        %SP(i).IP.mean_ssim_wn = mean(SP(i).IP.ssim_wn);
       [yPeaks,xPeaks] = findpeaks(SP(i).IP.ssim_wn, SP(i).wn, 'SortStr','descend');
       if (size(xPeaks) ~= 0)

           SP(i).IP.peaks_ssim_wn = [xPeaks(1) xPeaks(2) xPeaks(3)];
           %selection criteria
           val_temp = (yPeaks(1))^2 + ...
                       (yPeaks(2))^2 + ...
                       (yPeaks(3))^2;
       else
            SP(i).IP.peaks_ssim_wn = [0 0 0];
            val_temp = 0;
       end
        
       if(val_max < val_temp)
            val_max = val_temp;
            best_indice = i;
       end
        
    end
    best_ratio_temp = SP(best_indice).ratio_window;

    if(best_ratio_temp ~= 0.05)
        ratio_init = best_ratio_temp - 0.05;
    else
        ratio_init = 0.05;
    end
    if (best_ratio_temp ~= 1)
        ratio_fin = best_ratio_temp + 0.05;
    else
        ratio_fin = 1;
    end

    step = 0.01;
    ratios = ratio_init:step:ratio_fin;
    for i = 1:numel(ratios)
        SP(i) = SPtemp.copy(); %copying the data internally is much faster than loading all the files every time
        SP(i).ratio_window = ratios(i);
        SP(i) = SP(i).window_overlap_to_test(SP(i).tukey_window_param,SP(i).pourc_pulse_width);
        SP(i) = SP(i).Tnorm_and_center_data(1,0);
        SP(i) = SP(i).stitch_time_axis_T_with_interp(SP(i).interp_method);
        SP(i) = SP(i).pick_fourier_window(SP(i).window2_name); 
        SP(i) = SP(i).FT(SP(i).data_stitched.t_stitched, permute(SP(i).data_stitched.data_R,[3 1 2]).*repmat(SP(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
        SP(i) = SP(i).make_raman_spectrum();
    end
    
    best_indice = 1; %wanted interaction 
    val_temp = 0; 
    val_max = 0;
    for i = 1:numel(ratios)
        SP(i).IP = Image_Processing();

        temp2=squeeze(mean((SP(i).data_processed(1).data_T),[1]));
        temp2=(temp2-min(temp2(:)))./max(temp2(:));
        SP(i).IP.mat_ref = temp2;

        for j = 1:size(SP(i).hyperspectralRamanImageComplex,1)/2
            temp=squeeze(abs(SP(i).hyperspectralRamanImageComplex(j,:,:)));
            temp=(temp-min(temp(:)))./max(temp(:));
            SP(i).IP.mat_img_wn{j} = temp;
            SP(i).IP.ssim_wn(j) = ssim(temp,temp2);
        end
        %SP(i).IP.mean_ssim_wn = mean(SP(i).IP.ssim_wn);
       [yPeaks,xPeaks] = findpeaks(SP(i).IP.ssim_wn, SP(i).wn, 'SortStr','descend');
       if (size(xPeaks) ~= 0)

           SP(i).IP.peaks_ssim_wn = [xPeaks(1) xPeaks(2) xPeaks(3)];
           %selection criteria
           val_temp = (yPeaks(1))^2 + ...
                       (yPeaks(2))^2 + ...
                       (yPeaks(3))^2;
       else
            SP(i).IP.peaks_ssim_wn = [0 0 0];
            val_temp = 0;
       end
        
       if(val_max < val_temp)
            val_max = val_temp;
            best_indice = i;
       end
        
    end
    best_ratio_temp = SP(best_indice).ratio_window;
    fprintf('The best ratio window(%s) by ssim is %s.\n',SP(best_indice).window2_name,string(best_ratio_temp))
    best_ratio = best_ratio_temp;

end

