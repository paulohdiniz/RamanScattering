function get_best_filter(SP)
            %At least 9 min is needed for the complete report
            matrix_img_ref = SP.IP.mat_ref;
            cell_with_imgs_wn = SP.IP.mat_img_wn(SP.pixels_plot([1 2 3 4 5]));

            %TXT WITH RESULTS
            fileID = fopen('results_best_filter_carre.txt','w');
            fprintf(fileID, 'Analysis date: %s \n', datetime("now"));
            fprintf(fileID, 'Intern: Paulo Henrique DINIZ FERNANDES\n');
            fprintf(fileID, 'Tukey Window parameter: %f \n', SP.tukey_window_param);
            fprintf(fileID, 'Deadtime: %f \n', SP.deadtime);
            fprintf(fileID, 'Window type: %s \n', SP.window2_name);
            fprintf(fileID, 'Ratio Window: %f \n', SP.ratio_window);
            fprintf(fileID, 'Peaks in order: %d, %d, %d, %d, %d \n', SP.pixels_plot(1), SP.pixels_plot(2), SP.pixels_plot(3), SP.pixels_plot(4), SP.pixels_plot(5));
            fprintf(fileID, 'SSIM of the best images (peaks) without filter: %s, %s, %s, %s, %s \n', num2str(SP.IP.peaks_ssim(1)), num2str(SP.IP.peaks_ssim(2)), num2str(SP.IP.peaks_ssim(3)), num2str(SP.IP.peaks_ssim(4)), num2str(SP.IP.peaks_ssim(5)));
            fprintf(fileID, '-------------------------------------------\n');

            %AVG FILTER
            best_ssim_temp = -1; 
            best_i = 1;
            best_k = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i = 1:10
                    img_filtered_temp = SP.IP.average_filter(cell_with_imgs_wn{k}, [i i]);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_i = i;
                        best_k = k;
                    end
                end
            end
            fprintf(fileID, 'AVG_FILTER: The best ssim happens to image %d and is worth %f for the size filter [%d %d]. \n', best_k, best_ssim_temp, best_i, best_i);

            %MEDFILTER2
            best_ssim_temp = -1;
            best_i = 1;
            best_k = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=1:10
                    img_filtered_temp = SP.IP.medfilt2_filter(cell_with_imgs_wn{k}, i, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_i = i;
                        best_k = k;
                    end
                end
            end
            fprintf(fileID, 'MED2_FILTER: The best ssim happens to image %d and is worth %f for the size filter [%d %d]. \n', best_k, best_ssim_temp, best_i, best_i);

            %WIENER2
            best_ssim_temp = -1;
            best_i = 1;
            best_k = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=1:10
                    img_filtered_temp = SP.IP.wiener2_filter(cell_with_imgs_wn{k}, i, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_i = i;
                        best_k = k;
                    end
                end
            end
            fprintf(fileID, 'WIENER2_FILTER: The best ssim happens to image %d and is worth %f for the size filter [%d %d]. \n', best_k, best_ssim_temp, best_i, best_i);

            %GAUSSIEN1
            best_ssim_temp = -1;
            best_i = 1;
            best_k = 1;
            best_sigma = 0;
            for k=1:numel(cell_with_imgs_wn)
                for i=1:10
                    for sigma=0.1:0.05:5
                        img_filtered_temp = SP.IP.gaussian1_filter(cell_with_imgs_wn{k}, [i i], sigma);
                        ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                        if (ssim_temp > best_ssim_temp)
                            best_ssim_temp = ssim_temp;
                            best_i = i;
                            best_k = k;
                            best_sigma = sigma;
                        end
                    end
                end
            end
            fprintf(fileID, 'GAUSSIEN1_FILTER: The best ssim happens to image %d and is worth %f for the size filter [%d %d] with sigma %f. \n', best_k, best_ssim_temp, best_i, best_i, best_sigma);

            %GAUSSIEN2
            best_ssim_temp = -1;
            best_i = 1;
            best_k = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=0.1:0.05:5
                    img_filtered_temp = SP.IP.gaussian2_filter(cell_with_imgs_wn{k}, i, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_i = i;
                        best_k = k;
                    end
                end
            end
            fprintf(fileID, 'GAUSSIEN2_FILTER: The best ssim happens to image %d and is worth %f for the sigma [%f %f]. \n', best_k, best_ssim_temp, best_i, best_i);

            %MODE FILTER
            best_ssim_temp = -1;
            best_i = 1;
            best_k = 1;
            for k = 1:numel(cell_with_imgs_wn)
                for i = 1:2:11 %vector of positive odd integers
                    img_filtered_temp = SP.IP.mode_filter(cell_with_imgs_wn{k}, [i i]);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_i = i;
                        best_k = k;
                    end
                end
            end
            fprintf(fileID, 'MODE_FILTER: The best ssim happens to image %d and is worth %f for the size filter [%d %d]. \n', best_k, best_ssim_temp, best_i, best_i);
                 
            %DISK FILTER
            best_ssim_temp = -1;
            best_k = 1;
            best_i = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=2:11
                    img_filtered_temp = SP.IP.disk_filter(cell_with_imgs_wn{k}, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_k = k;
                        best_i = i;
                    end
                end
            end
            fprintf(fileID, 'DISK_FILTER: The best ssim happens to image %d and is worth %f with radius %d. \n', best_k, best_ssim_temp, best_i);
             
            %ORDFILT2 FILTER
            best_ssim_temp = -1;
            best_k = 1;
            best_i = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=2:9
                    img_filtered_temp = SP.IP.ordfilt2_filter(cell_with_imgs_wn{k}, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_k = k;
                        best_i = i;
                    end
                end
            end
            fprintf(fileID, 'ORDFILT2_FILTER: The best ssim happens to image %d and is worth %f with order %d. \n', best_k, best_ssim_temp, best_i);


            % imguided filter
            best_ssim_temp = -1;
            best_i = 1;
            best_k = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=1:10
                    img_filtered_temp = SP.IP.imguided_filter(cell_with_imgs_wn{k}, i, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_i = i;
                        best_k = k;
                    end
                end
            end
            fprintf(fileID, 'IMGUIDED_FILTER: The best ssim happens to image %d and is worth %f for the size filter [%d %d]. \n', best_k, best_ssim_temp, best_i, best_i);


            %LAPLACIEN FILTER
            best_ssim_temp = -1;
            best_k = 1;
            best_alpha = 1;
            for k=1:numel(cell_with_imgs_wn)
                for alpha=0:0.05:1
                    img_filtered_temp = SP.IP.laplacian_filter(cell_with_imgs_wn{k}, alpha);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_k = k;
                        best_alpha = alpha;
                    end
                end
            end
            fprintf(fileID, 'LAPLACIEN_FILTER: The best ssim happens to image %d and is worth %f with alpha %d. \n', best_k, best_ssim_temp, best_alpha);
             
            % imagesc(cell_with_imgs_wn{1})
            % imagesc(SP.IP.average_filter(cell_with_imgs_wn{1}, [5 7]))
            % a = SP.IP.average_filter(cell_with_imgs_wn{1}, [50 51]);
            % b = SP.IP.average_filter(cell_with_imgs_wn{1}, [50 50]);
            %selection criteria

            fclose(fileID);
            
            % val_temp = (yPeaks(1))^2 + ...
            %            (yPeaks(2))^2 + ...
            %            (yPeaks(3))^2 + ...
            %            (yPeaks(4))^2 + ...
            %             (yPeaks(5))^2;
  
end

