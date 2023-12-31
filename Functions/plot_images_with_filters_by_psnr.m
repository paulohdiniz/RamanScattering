function plot_images_with_filters_by_psnr(SP)
    %In the TESE SID: all the hyperspectral images in this section are 
    % smoothened with a Gaussian filter of standard deviation of 0.55.
    [t_values_for_stop, t_subsections, y_subsections, t, ampli, freq_of_peaks] = get_ideal_times_by_spectrogram(SP);
    
    for i=1:numel(t_values_for_stop)
        SPs(i) = SP.copy(); %copying the data internally is much faster than loading all the files every time
        SPs(i).ratio_window = t_values_for_stop(i)/max(SPs(i).data_stitched.t_stitched);
        SPs(i) = SPs(i).window_overlap_to_test(SPs(i).tukey_window_param,SPs(i).percent_FWHM);
        SPs(i) = SPs(i).Tnorm_and_center_data(1,0);
        SPs(i) = SPs(i).stitch_time_axis_T_with_interp(SPs(i).interp_method);
        SPs(i) = SPs(i).pick_fourier_window(SPs(i).window2_name); 
        SPs(i) = SPs(i).FT(SPs(i).data_stitched.t_stitched, permute(SPs(i).data_stitched.data_R,[3 1 2]).*repmat(SPs(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
        SPs(i) = SPs(i).make_raman_spectrum();
        SPs(i) = SPs(i).points_to_plot_by_frequency();   
    end

    matrix_img_ref = SP.IP.mat_ref;
    for i = 1:numel(SPs)
        peaks_imgs(i).img = squeeze(abs(SPs(i).hyperspectralRamanImageComplex(SP.pixels_plot(i),:,:)));
        peaks_imgs(i).img = peaks_imgs(i).img./max(peaks_imgs(i).img(:));
        peaks_imgs(i).piqe = piqe(peaks_imgs(i).img);
        peaks_imgs(i).niqe = niqe(peaks_imgs(i).img);
        peaks_imgs(i).brisque = brisque(peaks_imgs(i).img);
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
    filter(1).name = sprintf("AVG: [%d %d]", best_i, best_i);
    filter(1).psnr = best_psnr_temp;
    filter(1).img_t = best_img;
    for i=1:3
        filter(1).img_hyp(i).img = SP.IP.average_filter(peaks_imgs(i).img, [best_i,best_i]);
        filter(1).img_hyp(i).img = filter(1).img_hyp(i).img./max(filter(1).img_hyp(i).img(:));
        filter(1).img_hyp(i).piqe = piqe(filter(1).img_hyp(i).img);
        filter(1).img_hyp(i).niqe = niqe(filter(1).img_hyp(i).img);
        filter(1).img_hyp(i).brisque = brisque(filter(1).img_hyp(i).img);
        filter(1).img_hyp(i).psnr = psnr(filter(1).img_hyp(i).img, peaks_imgs(i).img);
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
    filter(2).name = sprintf("MED: [%d %d]", best_i, best_i);
    filter(2).psnr = best_psnr_temp;
    filter(2).img_t = best_img;
    for i=1:3
        filter(2).img_hyp(i).img = SP.IP.medfilt2_filter(peaks_imgs(i).img, best_i,best_i);
        filter(2).img_hyp(i).img = filter(2).img_hyp(i).img./max(filter(2).img_hyp(i).img(:));
        filter(2).img_hyp(i).piqe = piqe(filter(2).img_hyp(i).img);
        filter(2).img_hyp(i).niqe = niqe(filter(2).img_hyp(i).img);
        filter(2).img_hyp(i).brisque = brisque(filter(2).img_hyp(i).img);
        filter(2).img_hyp(i).psnr = psnr(filter(2).img_hyp(i).img, peaks_imgs(i).img);
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
    filter(3).name = sprintf("Wiener: [%d %d]", best_i, best_i);
    filter(3).psnr = best_psnr_temp;
    filter(3).img_t = best_img;
    for i=1:3
        filter(3).img_hyp(i).img = SP.IP.wiener2_filter(peaks_imgs(i).img, best_i,best_i);
        filter(3).img_hyp(i).img = filter(3).img_hyp(i).img./max(filter(3).img_hyp(i).img(:));
        filter(3).img_hyp(i).piqe = piqe(filter(3).img_hyp(i).img);
        filter(3).img_hyp(i).niqe = niqe(filter(3).img_hyp(i).img);
        filter(3).img_hyp(i).brisque = brisque(filter(3).img_hyp(i).img);
        filter(3).img_hyp(i).psnr = psnr(filter(3).img_hyp(i).img, peaks_imgs(i).img);
    end

    %GAUSSIEN1
    best_psnr_temp = -1;
    best_i = 1;
    best_img =[];
    best_sigma = 0.00;
    for i=2:10
        for sigma=0.5:0.05:5
            img_filtered_temp = SP.IP.gaussian1_filter(matrix_img_ref, [i i], sigma);
            psnr_temp = psnr(img_filtered_temp,matrix_img_ref);
             if (psnr_temp > best_psnr_temp)
                best_psnr_temp = psnr_temp;
                best_i = i;
                best_img = img_filtered_temp;
                best_sigma = sigma;
            end
        end
    end
    filter(4).name = sprintf("Gaussien1: [%d %d]; sigma = %.2f", best_i, best_i, best_sigma);
    filter(4).psnr = best_psnr_temp;
    filter(4).img_t = best_img;
    for i=1:3
        filter(4).img_hyp(i).img = SP.IP.gaussian1_filter(peaks_imgs(i).img, [best_i best_i], best_sigma);
        filter(4).img_hyp(i).img = filter(4).img_hyp(i).img./max(filter(4).img_hyp(i).img(:));
        filter(4).img_hyp(i).piqe = piqe(filter(4).img_hyp(i).img);
        filter(4).img_hyp(i).niqe = niqe(filter(4).img_hyp(i).img);
        filter(4).img_hyp(i).brisque = brisque(filter(4).img_hyp(i).img);
        filter(4).img_hyp(i).psnr = psnr(filter(4).img_hyp(i).img, peaks_imgs(i).img);
    end

    %DISK FILTER
    best_psnr_temp = -1;
    best_i = 1;
    best_img =[];
    for i=2:11
        img_filtered_temp = SP.IP.disk_filter(matrix_img_ref, i);
        psnr_temp = psnr(img_filtered_temp,matrix_img_ref);
         if (psnr_temp > best_psnr_temp)
            best_psnr_temp = psnr_temp;
            best_i = i;
            best_img = img_filtered_temp;
        end
    end
    filter(5).name = sprintf("Disk: %d radius", best_i);
    filter(5).psnr = best_psnr_temp;
    filter(5).img_t = best_img;
    for i=1:3
        filter(5).img_hyp(i).img = SP.IP.disk_filter(peaks_imgs(i).img,best_i);
        filter(5).img_hyp(i).img = filter(5).img_hyp(i).img./max(filter(5).img_hyp(i).img(:));
        filter(5).img_hyp(i).piqe = piqe(filter(5).img_hyp(i).img);
        filter(5).img_hyp(i).niqe = niqe(filter(5).img_hyp(i).img);
        filter(5).img_hyp(i).brisque = brisque(filter(5).img_hyp(i).img);
        filter(5).img_hyp(i).psnr = psnr(filter(5).img_hyp(i).img, peaks_imgs(i).img);
    end

    %ORDFILT2 FILTER
    best_psnr_temp = -1;
    best_i = 1;
    best_img =[];
    for i=2:9
        img_filtered_temp = SP.IP.ordfilt2_filter(matrix_img_ref, i);
        psnr_temp = psnr(img_filtered_temp,matrix_img_ref);
         if (psnr_temp > best_psnr_temp)
            best_psnr_temp = psnr_temp;
            best_i = i;
            best_img = img_filtered_temp;
        end
    end
    filter(6).name = sprintf("ORD: %d", best_i);
    filter(6).psnr = best_psnr_temp;
    filter(6).img_t = best_img;
    for i=1:3
        filter(6).img_hyp(i).img = SP.IP.ordfilt2_filter(peaks_imgs(i).img,best_i);
        filter(6).img_hyp(i).img = filter(6).img_hyp(i).img./max(filter(6).img_hyp(i).img(:));
        filter(6).img_hyp(i).piqe = piqe(filter(6).img_hyp(i).img);
        filter(6).img_hyp(i).niqe = niqe(filter(6).img_hyp(i).img);
        filter(6).img_hyp(i).brisque = brisque(filter(6).img_hyp(i).img);
        filter(6).img_hyp(i).psnr = psnr(filter(6).img_hyp(i).img, peaks_imgs(i).img);
    end

    % imguided filter
    best_psnr_temp = -1;
    best_i = 1;
    best_img =[];
    for i=2:10
        img_filtered_temp = SP.IP.imguided_filter(matrix_img_ref, i, i);
        psnr_temp = psnr(img_filtered_temp,matrix_img_ref);
         if (psnr_temp > best_psnr_temp)
            best_psnr_temp = psnr_temp;
            best_i = i;
            best_img = img_filtered_temp;
        end
    end
    filter(7).name = sprintf("Imguided: [%d %d]", best_i,best_i);
    filter(7).psnr = best_psnr_temp;
    filter(7).img_t = best_img;
    for i=1:3
        filter(7).img_hyp(i).img = SP.IP.imguided_filter(peaks_imgs(i).img,best_i, best_i);
        filter(7).img_hyp(i).img = filter(7).img_hyp(i).img./max(filter(7).img_hyp(i).img(:));
        filter(7).img_hyp(i).piqe = piqe(filter(7).img_hyp(i).img);
        filter(7).img_hyp(i).niqe = niqe(filter(7).img_hyp(i).img);
        filter(7).img_hyp(i).brisque = brisque(filter(7).img_hyp(i).img);
        filter(7).img_hyp(i).psnr = psnr(filter(7).img_hyp(i).img, peaks_imgs(i).img);
    end
 
    name_of_figure = append('Exp : ', string(SP.xp_number));
    h1 = figure('Position', [50 100 1300 800],'Name', name_of_figure);
    
    for i=1:numel(filter)
        psnr_vec = [filter(:).psnr];
        name_vec = {filter(:).name};
        img_hyp_vec = {filter(:).img_hyp};
        img_t_vec = {filter(:).img_t};

        [~, sortedIndices] = sort(psnr_vec, 'descend');

        topThreeIndices = sortedIndices(1:3);   

        topThree_img_t = img_t_vec(topThreeIndices);
        topThree_psnr = psnr_vec(topThreeIndices);
        topThree_name = name_vec(topThreeIndices);
        topThree_img_hyp = img_hyp_vec(topThreeIndices);
    end
        
        %% PLOTS
        subplot(5,9,[1,2])
        imagesc(matrix_img_ref);
        axis image
        xlabel('pixels','interpreter','tex');
        ylabel({[append('niqe: ', num2str(niqe(matrix_img_ref),'%.2f'))] ...
            [append('brisque: ', num2str(brisque(matrix_img_ref),'%.2f'))] ...
            [append('piqe: ', num2str(piqe(matrix_img_ref),'%.2f'))]}, ...
        'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','normal', ...
        'fontsize',9, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom');
        yticks([])
        xticks([0 10 20 30 40 50])
         
        title('Transmission Image','fontsize',12)
        text(-0.1, 1.1,'','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)

        subplot(5,9,[4,5])
        imagesc(topThree_img_t{1});
        axis image
        title(topThree_name{1},'fontsize',10);
        subtitle(append('psnr: ',num2str(topThree_psnr(1),'%.2f')))
        ylabel({[append('niqe: ', num2str(niqe(topThree_img_t{1}),'%.2f'))] ...
            [append('brisque: ', num2str(brisque(topThree_img_t{1}),'%.2f'))] ...
            [append('piqe: ', num2str(piqe(topThree_img_t{1}),'%.2f'))] ...
            [append('{\Delta}P: ',num2str(100*(piqe(topThree_img_t{1})- piqe(matrix_img_ref))/piqe(matrix_img_ref),'%.2f'),'%')]}, ...
        'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','normal', ...
        'fontsize',9, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom',...
        'interpreter','tex');
        xticks([]);
        yticks([]);
        colormap('hot')

        subplot(5,9,[6,7])
        imagesc(topThree_img_t{2});
        axis image
        title(topThree_name{2},'fontsize',10);
        subtitle(append('psnr: ',num2str(topThree_psnr(2),'%.2f')))
        ylabel({[append('niqe: ', num2str(niqe(topThree_img_t{2}),'%.2f'))] ...
            [append('brisque: ', num2str(brisque(topThree_img_t{2}),'%.2f'))] ...
            [append('piqe: ', num2str(piqe(topThree_img_t{2}),'%.2f'))] ...
            [append('{\Delta}P: ',num2str(100*(piqe(topThree_img_t{2})- piqe(matrix_img_ref))/piqe(matrix_img_ref),'%.2f'),'%')]}, ...
        'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','normal', ...
        'fontsize',9, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom',...
        'interpreter','tex');
        xticks([]);
        yticks([]);
        colormap('hot')

        subplot(5,9,[8,9])
        imagesc(topThree_img_t{3});
        axis image
        title(topThree_name{3},'fontsize',10);
        subtitle(append('psnr: ',num2str(topThree_psnr(3),'%.2f')))
        ylabel({[append('niqe: ', num2str(niqe(topThree_img_t{3}),'%.2f'))] ...
            [append('brisque: ', num2str(brisque(topThree_img_t{3}),'%.2f'))] ...
            [append('piqe: ', num2str(piqe(topThree_img_t{3}),'%.2f'))] ...
            [append('{\Delta}P: ',num2str(100*(piqe(topThree_img_t{3})- piqe(matrix_img_ref))/piqe(matrix_img_ref),'%.2f'),'%')]}, ...
        'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','normal', ...
        'fontsize',9, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom',...
        'interpreter','tex');
        xticks([]);
        yticks([]);        
        colormap('hot')

        subplot(5,9, [19, 20])
        imagesc(peaks_imgs(1).img);
        axis image
        title(append(num2str(SP.wns_plot(1), '%.2f'), ' cm^{-1}'),'fontsize',10);
        ylabel({[append('niqe: ', num2str(peaks_imgs(1).niqe,'%.2f'))] ...
            [append('brisque: ', num2str(peaks_imgs(1).brisque,'%.2f'))] ...
            [append('piqe: ', num2str(peaks_imgs(1).piqe,'%.2f'))]}, 'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','bold', ...
        'fontsize',9, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom',...
        'interpreter','tex');
        xticks([]);
        yticks([]);
        colormap('hot')

        subplot(5,9, [28, 29])
        imagesc(peaks_imgs(2).img);
        axis image
        title(append(num2str(SP.wns_plot(2), '%.2f'), ' cm^{-1}'),'fontsize',10);
        ylabel({[append('niqe: ', num2str(peaks_imgs(2).niqe,'%.2f'))] ...
            [append('brisque: ', num2str(peaks_imgs(2).brisque,'%.2f'))] ...
            [append('piqe: ', num2str(peaks_imgs(2).piqe,'%.2f'))]}, 'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','bold', ...
        'fontsize',9, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom',...
        'interpreter','tex');
        xticks([]);
        yticks([]);
        colormap('hot')

        subplot(5,9, [37, 38])
        imagesc(peaks_imgs(3).img);
        axis image
        title(append(num2str(SP.wns_plot(3), '%.2f'), ' cm^{-1}'),'fontsize',10);
        ylabel({[append('niqe: ', num2str(peaks_imgs(3).niqe,'%.2f'))] ...
            [append('brisque: ', num2str(peaks_imgs(3).brisque,'%.2f'))] ...
            [append('piqe: ', num2str(peaks_imgs(3).piqe,'%.2f'))]}, 'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','bold', ...
        'fontsize',9, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom',...
        'interpreter','tex');
        xticks([]);
        yticks([]);
        colormap('hot')

        for k = [3, 4, 5]
    
            subplot(5,9,[9*k-5,9*k-4])
            imagesc(topThree_img_hyp{1}(k-2).img);
            axis image
            subtitle(append('psnr: ',num2str(psnr(topThree_img_hyp{1}(k-2).img,peaks_imgs(k-2).img),'%.2f')))
            ylabel({[append('niqe: ', num2str(niqe(topThree_img_hyp{1}(k-2).img),'%.2f'))] ...
                [append('brisque: ', num2str(brisque(topThree_img_hyp{1}(k-2).img),'%.2f'))] ...
                [append('piqe: ', num2str(piqe(topThree_img_hyp{1}(k-2).img),'%.2f'))] ...
                [append('{\Delta}P: ',num2str(100*(piqe(topThree_img_hyp{1}(k-2).img) - piqe(peaks_imgs(k-2).img))/piqe(peaks_imgs(k-2).img),'%.2f'),'%')]}, ...
                'Rotation',0, ...
                'interpreter','none', ...
                'fontweight','normal', ...
                'fontsize',9, ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','bottom',...
                'interpreter','tex');
            xticks([]);
            yticks([]);
            colormap('hot')

            subplot(5,9,[9*k-3,9*k-2])
            imagesc(topThree_img_hyp{2}(k-2).img);
            axis image
            subtitle(append('psnr: ',num2str(psnr(topThree_img_hyp{2}(k-2).img,peaks_imgs(k-2).img),'%.2f')))
            ylabel({[append('niqe: ', num2str(niqe(topThree_img_hyp{2}(k-2).img),'%.2f'))] ...
                [append('brisque: ', num2str(brisque(topThree_img_hyp{2}(k-2).img),'%.2f'))] ...
                [append('piqe: ', num2str(piqe(topThree_img_hyp{2}(k-2).img),'%.2f'))] ...
                [append('{\Delta}P: ',num2str(100*(piqe(topThree_img_hyp{2}(k-2).img) - piqe(peaks_imgs(k-2).img))/piqe(peaks_imgs(k-2).img),'%.2f'),'%')]}, ...
                'Rotation',0, ...
                'interpreter','none', ...
                'fontweight','normal', ...
                'fontsize',9, ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','bottom',...
                'interpreter','tex');
            xticks([]);
            yticks([]);
            colormap('hot')

            subplot(5,9,[9*k-1,9*k])
            imagesc(topThree_img_hyp{3}(k-2).img);
            axis image
            subtitle(append('psnr: ',num2str(psnr(topThree_img_hyp{3}(k-2).img,peaks_imgs(k-2).img),'%.2f')))
            ylabel({[append('niqe: ', num2str(niqe(topThree_img_hyp{3}(k-2).img),'%.2f'))] ...
                [append('brisque: ', num2str(brisque(topThree_img_hyp{3}(k-2).img),'%.2f'))] ...
                [append('piqe: ', num2str(piqe(topThree_img_hyp{3}(k-2).img),'%.2f'))] ...
                [append('{\Delta}P: ',num2str(100*(piqe(topThree_img_hyp{3}(k-2).img) - piqe(peaks_imgs(k-2).img))/piqe(peaks_imgs(k-2).img),'%.2f'),'%')]}, ...
                'Rotation',0, ...
                'interpreter','none', ...
                'fontweight','normal', ...
                'fontsize',9, ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','bottom',...
                'interpreter','tex');
            xticks([]);
            yticks([]);
            colormap('hot')
        end

end
