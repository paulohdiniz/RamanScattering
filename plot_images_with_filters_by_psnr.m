function plot_images_with_filters_by_psnr(SP)

    [t_values_for_stop, t_subsections, y_subsections, t, ampli, freq_of_peaks] = get_ideal_times_by_spectrogram(SP);
    
    for i=1:numel(t_values_for_stop)
        SPs(i) = SP.copy(); %copying the data internally is much faster than loading all the files every time
        SPs(i).ratio_window = t_values_for_stop(i)/max(SPs(i).data_stitched.t_stitched);
        SPs(i) = SPs(i).window_overlap_to_test(SPs(i).tukey_window_param,SPs(i).deadtime);
        SPs(i) = SPs(i).Tnorm_and_center_data(1,0,SPs(i).deadtime);
        SPs(i) = SPs(i).stitch_time_axis_T_with_interp(SPs(i).interp_method);
        SPs(i) = SPs(i).pick_fourier_window(SPs(i).window2_name); 
        SPs(i) = SPs(i).FT(SPs(i).data_stitched.t_stitched, permute(SPs(i).data_stitched.data_R,[3 1 2]).*repmat(SPs(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
        SPs(i) = SPs(i).make_raman_spectrum();
        SPs(i) = SPs(i).points_to_plot_by_frequency();   
    end

    matrix_img_ref = SP.IP.mat_ref;
    for i = 1:numel(SPs)
        cell_with_imgs{i} = squeeze(abs(SPs(i).hyperspectralRamanImageComplex(SP.pixels_plot(i),:,:)));
    end

    %AVG FILTER
    best_psnr_temp = -1; 
    best_i = 1;
    best_img =[];
    for i = 2:10
        img_filtered_temp = SP.IP.average_filter(matrix_img_ref, [i i]);
        psnr_temp = psnr(img_filtered_temp,matrix_img_ref);
        if (psnr_temp > best_psnr_temp)
            best_psnr_temp = psnr_temp;
            best_i = i;
            best_img = img_filtered_temp;
        end
    end
    filter(1).name = "AVG";
    filter(1).psnr = best_psnr_temp;
    filter(1).i = best_i;
    filter(1).img_t = best_img;
    for i=1:3
        filter(1).img_hyp(i).img = SP.IP.average_filter(cell_with_imgs{i}, [best_i,best_i]);
    end

    %MEDFILTER2
    best_psnr_temp = -1;
    best_i = 1;
    best_img =[];
    for i=2:10
        img_filtered_temp = SP.IP.medfilt2_filter(matrix_img_ref, i, i);
        psnr_temp = psnr(img_filtered_temp,matrix_img_ref);
        if (psnr_temp > best_psnr_temp)
            best_psnr_temp = psnr_temp;
            best_i = i;
            best_img = img_filtered_temp;
        end
    end
    filter(2).name = "MED";
    filter(2).psnr = best_psnr_temp;
    filter(2).i = best_i;
    filter(2).img_t = best_img;
    for i=1:3
        filter(2).img_hyp(i).img = SP.IP.medfilt2_filter(cell_with_imgs{i}, best_i,best_i);
    end

    %WIENER2
    best_psnr_temp = -1;
    best_i = 1;
    best_img =[];
    for i=2:10
        img_filtered_temp = SP.IP.wiener2_filter(matrix_img_ref, i, i);
        psnr_temp = psnr(img_filtered_temp,matrix_img_ref);
        if (psnr_temp > best_psnr_temp)
            best_psnr_temp = psnr_temp;
            best_i = i;
            best_img = img_filtered_temp;
        end
    end
    filter(3).name = "WIENER";
    filter(3).psnr = best_psnr_temp;
    filter(3).i = best_i;
    filter(3).img_t = best_img;
    for i=1:3
        filter(3).img_hyp(i).img = SP.IP.wiener2_filter(cell_with_imgs{i}, best_i,best_i);
    end

    name_of_figure = append('Exp : ', string(SP.xp_number));
    h1 = figure('Position', [50 100 1300 800],'Name', name_of_figure);
    for i=1:numel(filter)
        img_t_vec = {filter(:).img_t};
        psnr_vec = [filter(:).psnr];
        name_vec = {filter(:).name};
        param_vec = {filter(:).i};
        img_hyp_vec = {filter(:).img_hyp};

        [~, sortedIndices] = sort(psnr_vec, 'descend');

        topThreeIndices = sortedIndices(1:3);   

        topThree_img_t = img_t_vec(topThreeIndices);
        topThree_psnr = psnr_vec(topThreeIndices);
        topThree_name = name_vec(topThreeIndices);
        topThree_param = param_vec(topThreeIndices);
        topThree_img_hyp = img_hyp_vec(topThreeIndices);
    end
        
        %% PLOTS
        subplot(5,9,[1,2])
        imagesc(squeeze(mean((SP.data_processed(1).data_T),[1])));
        axis image
        title('DC PD signal','fontsize',14); 
        xlabel('pixels','interpreter','tex');
        ylabel('pixels');
        colorbar;
        title('Transmission Image','fontsize',14)
        text(-0.1, 1.1,'','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)

        subplot(5,9,[4,5])
        imagesc(topThree_img_t{1});
        axis image
        xlabel('pixels','fontsize',6);
        ylabel('pixels','fontsize',6);
        colorbar;
        % title({
        %     [topThree_name{1} ' ' '[' num2str(topThree_param{1}) ']'] 
        %     [ 'SSIM: ' num2str(SPs(i).IP.ssim_wn(SP.pixels_plot(i)))]
        %     },'fontsize',8);
        colormap('hot')

        subplot(5,9,[6,7])
        imagesc(topThree_img_t{2});
        axis image
        xlabel('pixels','fontsize',6);
        ylabel('pixels','fontsize',6);
        colorbar;
        % title({
        %     [topThree_name{2} ' ' '[' num2str(topThree_param{2}) ']'] 
        %     [ 'SSIM: ' num2str(SPs(i).IP.ssim_wn(SP.pixels_plot(i)))]
        %     },'fontsize',8);
        colormap('hot')

        subplot(5,9,[8,9])
        imagesc(topThree_img_t{3});
        axis image
        xlabel('pixels','fontsize',6);
        ylabel('pixels','fontsize',6);
        colorbar;
        % title({
        %     [topThree_name{3} ' ' '[' num2str(topThree_param{3}) ']'] 
        %     [ 'SSIM: ' num2str(SPs(i).IP.ssim_wn(SP.pixels_plot(i)))]
        %     },'fontsize',8);
        colormap('hot')

        subplot(5,9, [19, 20])
        imagesc(cell_with_imgs{1});
        axis image
        xlabel('pixels','fontsize',6);
        ylabel('pixels','fontsize',6);
        colorbar;
        % title({
        %     [topThree_name{3} ' ' '[' num2str(topThree_param{3}) ']'] 
        %     [ 'SSIM: ' num2str(SPs(i).IP.ssim_wn(SP.pixels_plot(i)))]
        %     },'fontsize',8);
        colormap('hot')

        subplot(5,9, [28, 29])
        imagesc(cell_with_imgs{2});
        axis image
        xlabel('pixels','fontsize',6);
        ylabel('pixels','fontsize',6);
        colorbar;
        % title({
        %     [topThree_name{3} ' ' '[' num2str(topThree_param{3}) ']'] 
        %     [ 'SSIM: ' num2str(SPs(i).IP.ssim_wn(SP.pixels_plot(i)))]
        %     },'fontsize',8);
        colormap('hot')

        subplot(5,9, [37, 38])
        imagesc(cell_with_imgs{3});
        axis image
        xlabel('pixels','fontsize',6);
        ylabel('pixels','fontsize',6);
        colorbar;
        % title({
        %     [topThree_name{3} ' ' '[' num2str(topThree_param{3}) ']'] 
        %     [ 'SSIM: ' num2str(SPs(i).IP.ssim_wn(SP.pixels_plot(i)))]
        %     },'fontsize',8);
        colormap('hot')

        for k = [3, 4, 5]
    
            subplot(5,9,[9*k-5,9*k-4])
            imagesc(topThree_img_hyp{1}(k-2).img);
            axis image
            xlabel('pixels','fontsize',6);
            ylabel('pixels','fontsize',6);
            colorbar;
            % title({
            %     [topThree_name{3} ' ' '[' num2str(topThree_param{3}) ']'] 
            %     [ 'SSIM: ' num2str(SPs(i).IP.ssim_wn(SP.pixels_plot(i)))]
            %     },'fontsize',8);
            colormap('hot')

            subplot(5,9,[9*k-3,9*k-2])
            imagesc(topThree_img_hyp{2}(k-2).img);
            axis image
            xlabel('pixels','fontsize',6);
            ylabel('pixels','fontsize',6);
            colorbar;
            % title({
            %     [topThree_name{3} ' ' '[' num2str(topThree_param{3}) ']'] 
            %     [ 'SSIM: ' num2str(SPs(i).IP.ssim_wn(SP.pixels_plot(i)))]
            %     },'fontsize',8);
            colormap('hot')

            subplot(5,9,[9*k-1,9*k])
            imagesc(topThree_img_hyp{3}(k-2).img);
            axis image
            xlabel('pixels','fontsize',6);
            ylabel('pixels','fontsize',6);
            colorbar;
            % title({
            %     [topThree_name{3} ' ' '[' num2str(topThree_param{3}) ']'] 
            %     [ 'SSIM: ' num2str(SPs(i).IP.ssim_wn(SP.pixels_plot(i)))]
            %     },'fontsize',8);
            colormap('hot')
        end



end
