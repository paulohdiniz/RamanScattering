function get_txt_best_filter_by_ssim(SP)
            matrix_img_ref = SP.IP.mat_ref;
            cell_with_imgs_wn = SP.IP.mat_img_wn(SP.pixels_plot([1 2 3]));

            %TXT WITH RESULTS
            fileID = fopen(append('results_best_filter_byssim_exp',string(SP.xp_number), '.txt'),'w');
            fprintf(fileID, 'Analysis date: %s \n', datetime("now"));
            fprintf(fileID, 'Intern: Paulo Henrique DINIZ FERNANDES\n');
            fprintf(fileID, 'Folder path : %s \n', string(pwd));
            fprintf(fileID, 'Experiment number: %d \n', SP.xp_number);
            fprintf(fileID, 'Tukey Window parameter: %.4f \n', SP.tukey_window_param);
            fprintf(fileID, 'Deadtime: %.4f \n', SP.deadtime);
            fprintf(fileID, 'Window type: %s \n', SP.window2_name);
            fprintf(fileID, 'Ratio Window: %.4f \n', SP.ratio_window);
            fprintf(fileID, 'Peaks in order: %d, %d, %d, %d, %d \n', SP.pixels_plot(1), SP.pixels_plot(2), SP.pixels_plot(3));
            fprintf(fileID, 'SSIM of peaks without filter: %s, %s, %s \n', num2str(SP.IP.ssim_wn(SP.pixels_plot(1))), num2str(SP.IP.ssim_wn(SP.pixels_plot(2))), num2str(SP.IP.ssim_wn(SP.pixels_plot(3))));
            fprintf(fileID, 'BRISQUE of peaks without filter: %s, %s, %s \n', num2str(SP.IP.brisque_score_wn(SP.pixels_plot(1))), num2str(SP.IP.brisque_score_wn(SP.pixels_plot(2))), num2str(SP.IP.brisque_score_wn(SP.pixels_plot(3))));
            fprintf(fileID, 'NIQE of peaks without filter: %s, %s, %s \n', num2str(SP.IP.niqe_score_wn(SP.pixels_plot(1))), num2str(SP.IP.niqe_score_wn(SP.pixels_plot(2))), num2str(SP.IP.niqe_score_wn(SP.pixels_plot(3))));
            fprintf(fileID, 'PIQE of peaks without filter: %s, %s, %s \n', num2str(SP.IP.piqe_score_wn(SP.pixels_plot(1))), num2str(SP.IP.piqe_score_wn(SP.pixels_plot(2))), num2str(SP.IP.piqe_score_wn(SP.pixels_plot(3))));

            
            fprintf(fileID, '-------------------------------------------\n');

            %AVG FILTER
            best_ssim_temp = -1; 
            delta_ssim_temp = 0;
            best_delta_ssim_temp = 0;
            best_i = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i = 1:10
                    img_filtered_temp = SP.IP.average_filter(cell_with_imgs_wn{k}, [i i]);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > SP.IP.ssim_wn(SP.pixels_plot(k)))
                        delta_ssim_temp = (ssim_temp - SP.IP.ssim_wn(SP.pixels_plot(k)));
                        if (delta_ssim_temp > best_delta_ssim_temp)
                            best_ssim_temp = ssim_temp;
                            best_delta_ssim_temp = delta_ssim_temp;
                            best_i = i;
                        end
                    end
                end
                fprintf(fileID, 'AVG_FILTER: The best delta_ssim to image %d has the size filter [%d %d]. Before SSIM: %.4f. After SSIM: %.4f. Diff = %.2f%% \n', k, best_i, best_i, SP.IP.ssim_wn(SP.pixels_plot(k)), best_ssim_temp, best_delta_ssim_temp*100/SP.IP.ssim_wn(SP.pixels_plot(k)));
                best_ssim_temp = -1; 
                best_delta_ssim_temp = 0;
            end

            %MEDFILTER2
            best_ssim_temp = -1;
            delta_ssim_temp = 0;
            best_delta_ssim_temp = 0;
            best_i = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=1:10
                    img_filtered_temp = SP.IP.medfilt2_filter(cell_with_imgs_wn{k}, i, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > SP.IP.ssim_wn(SP.pixels_plot(k)))
                        delta_ssim_temp = (ssim_temp - SP.IP.ssim_wn(SP.pixels_plot(k)));
                        if (delta_ssim_temp > best_delta_ssim_temp)
                            best_ssim_temp = ssim_temp;
                            best_delta_ssim_temp = delta_ssim_temp;
                            best_i = i;
                        end
                    end
                end
                fprintf(fileID, 'MED2_FILTER: The best delta_ssim to image %d has the size filter [%d %d]. Before SSIM: %.4f. After SSIM: %.4f. Diff = %.2f%% \n', k, best_i, best_i, SP.IP.ssim_wn(SP.pixels_plot(k)), best_ssim_temp, best_delta_ssim_temp*100/SP.IP.ssim_wn(SP.pixels_plot(k)));
                best_ssim_temp = -1; 
                best_delta_ssim_temp = 0;
            end

            %WIENER2
            best_ssim_temp = -1;
            delta_ssim_temp = 0;
            best_delta_ssim_temp = 0;
            best_i = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=1:10
                    img_filtered_temp = SP.IP.wiener2_filter(cell_with_imgs_wn{k}, i, i);
                     ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > SP.IP.ssim_wn(SP.pixels_plot(k)))
                        delta_ssim_temp = (ssim_temp - SP.IP.ssim_wn(SP.pixels_plot(k)));
                        if (delta_ssim_temp > best_delta_ssim_temp)
                            best_ssim_temp = ssim_temp;
                            best_delta_ssim_temp = delta_ssim_temp;
                            best_i = i;
                        end
                    end
                end
                fprintf(fileID, 'WIENER2_FILTER: The best delta_ssim to image %d has the size filter [%d %d]. Before SSIM: %.4f. After SSIM: %.4f. Diff = %.2f%% \n', k, best_i, best_i, SP.IP.ssim_wn(SP.pixels_plot(k)), best_ssim_temp, best_delta_ssim_temp*100/SP.IP.ssim_wn(SP.pixels_plot(k)));
                best_ssim_temp = -1; 
                best_delta_ssim_temp = 0;
            end

            %GAUSSIEN1
            best_ssim_temp = -1;
            delta_ssim_temp = 0;
            best_delta_ssim_temp = 0;
            best_i = 1;
            best_sigma = 0;
            for k=1:numel(cell_with_imgs_wn)
                for i=1:10
                    for sigma=0.1:0.05:5
                        img_filtered_temp = SP.IP.gaussian1_filter(cell_with_imgs_wn{k}, [i i], sigma);
                        ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                        if (ssim_temp > SP.IP.ssim_wn(SP.pixels_plot(k)))
                            delta_ssim_temp = (ssim_temp - SP.IP.ssim_wn(SP.pixels_plot(k)));
                            if (delta_ssim_temp > best_delta_ssim_temp)
                                best_ssim_temp = ssim_temp;
                                best_delta_ssim_temp = delta_ssim_temp;
                                best_sigma = sigma;
                                best_i = i;
                            end
                        end
                    end
                end
                fprintf(fileID, 'GAUSSIEN1_FILTER: The best delta_ssim to image %d has the size filter [%d %d] and sigma %.2f. Before SSIM: %.4f. After SSIM: %.4f. Diff = %.2f%% \n', k, best_i, best_i, best_sigma, SP.IP.ssim_wn(SP.pixels_plot(k)), best_ssim_temp, best_delta_ssim_temp*100/SP.IP.ssim_wn(SP.pixels_plot(k)));
                best_ssim_temp = -1; 
                best_delta_ssim_temp = 0;
            end

            %GAUSSIEN2
            best_ssim_temp = -1;
            delta_ssim_temp = 0;
            best_delta_ssim_temp = 0;
            best_i = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=0.1:0.05:5
                    img_filtered_temp = SP.IP.gaussian2_filter(cell_with_imgs_wn{k}, i, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > SP.IP.ssim_wn(SP.pixels_plot(k)))
                        delta_ssim_temp = (ssim_temp - SP.IP.ssim_wn(SP.pixels_plot(k)));
                        if (delta_ssim_temp > best_delta_ssim_temp)
                            best_ssim_temp = ssim_temp;
                            best_delta_ssim_temp = delta_ssim_temp;
                            best_i = i;
                        end
                    end
                end
                fprintf(fileID, 'GAUSSIEN2_FILTER: The best delta_ssim to image %d has the sigma [%.2f %.2f]. Before SSIM: %.4f. After SSIM: %.4f. Diff = %.2f%% \n',  k, best_i, best_i, SP.IP.ssim_wn(SP.pixels_plot(k)), best_ssim_temp, best_delta_ssim_temp*100/SP.IP.ssim_wn(SP.pixels_plot(k)));
                best_ssim_temp = -1; 
                best_delta_ssim_temp = 0;
            end

            %MODE FILTER
            best_ssim_temp = -1;
            delta_ssim_temp = 0;
            best_delta_ssim_temp = 0;
            best_i = 1;
            for k = 1:numel(cell_with_imgs_wn)
                for i = 1:2:11 %vector of positive odd integers
                    img_filtered_temp = SP.IP.mode_filter(cell_with_imgs_wn{k}, [i i]);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > SP.IP.ssim_wn(SP.pixels_plot(k)))
                        delta_ssim_temp = (ssim_temp - SP.IP.ssim_wn(SP.pixels_plot(k)));
                        if (delta_ssim_temp > best_delta_ssim_temp)
                            best_ssim_temp = ssim_temp;
                            best_delta_ssim_temp = delta_ssim_temp;
                            best_i = i;
                        end
                    end
                end
                fprintf(fileID, 'MODE_FILTER: The best delta_ssim to image %d has the size filter [%d %d]. Before SSIM: %.4f. After SSIM: %.4f. Diff = %.2f%% \n',  k, best_i, best_i, SP.IP.ssim_wn(SP.pixels_plot(k)), best_ssim_temp, best_delta_ssim_temp*100/SP.IP.ssim_wn(SP.pixels_plot(k)));
                best_ssim_temp = -1; 
                best_delta_ssim_temp = 0;
            end
       
            %DISK FILTER
            best_ssim_temp = -1;
            delta_ssim_temp = 0;
            best_delta_ssim_temp = 0;
            best_i = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=2:11
                    img_filtered_temp = SP.IP.disk_filter(cell_with_imgs_wn{k}, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > SP.IP.ssim_wn(SP.pixels_plot(k)))
                        delta_ssim_temp = (ssim_temp - SP.IP.ssim_wn(SP.pixels_plot(k)));
                        if (delta_ssim_temp > best_delta_ssim_temp)
                            best_ssim_temp = ssim_temp;
                            best_delta_ssim_temp = delta_ssim_temp;
                            best_i = i;
                        end
                    end
                end
                fprintf(fileID, 'DISK_FILTER: The best delta_ssim to image %d has radius %d. Before SSIM: %.4f. After SSIM: %.4f. Diff = %.2f%% \n', k, best_i, SP.IP.ssim_wn(SP.pixels_plot(k)), best_ssim_temp, best_delta_ssim_temp*100/SP.IP.ssim_wn(SP.pixels_plot(k)));
                best_ssim_temp = -1; 
                best_delta_ssim_temp = 0;
            end
     
            %ORDFILT2 FILTER
            best_ssim_temp = -1;
            delta_ssim_temp = 0;
            best_delta_ssim_temp = 0;
            best_i = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=2:9
                    img_filtered_temp = SP.IP.ordfilt2_filter(cell_with_imgs_wn{k}, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > SP.IP.ssim_wn(SP.pixels_plot(k)))
                        delta_ssim_temp = (ssim_temp - SP.IP.ssim_wn(SP.pixels_plot(k)));
                        if (delta_ssim_temp > best_delta_ssim_temp)
                            best_ssim_temp = ssim_temp;
                            best_delta_ssim_temp = delta_ssim_temp;
                            best_i = i;
                        end
                    end
                end
                fprintf(fileID, 'ORDFILT2_FILTER: The best delta_ssim to image %d has order %d. Before SSIM: %.4f. After SSIM: %.4f. Diff = %.2f%% \n', k, best_i, SP.IP.ssim_wn(SP.pixels_plot(k)), best_ssim_temp, best_delta_ssim_temp*100/SP.IP.ssim_wn(SP.pixels_plot(k)));
                best_ssim_temp = -1; 
                best_delta_ssim_temp = 0;
            end

            % imguided filter
            best_ssim_temp = -1;
            delta_ssim_temp = 0;
            best_delta_ssim_temp = 0;
            best_i = 1;
            for k=1:numel(cell_with_imgs_wn)
                for i=1:10
                    img_filtered_temp = SP.IP.imguided_filter(cell_with_imgs_wn{k}, i, i);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > SP.IP.ssim_wn(SP.pixels_plot(k)))
                        delta_ssim_temp = (ssim_temp - SP.IP.ssim_wn(SP.pixels_plot(k)));
                        if (delta_ssim_temp > best_delta_ssim_temp)
                            best_ssim_temp = ssim_temp;
                            best_delta_ssim_temp = delta_ssim_temp;
                            best_i = i;
                        end
                    end
                end
                fprintf(fileID, 'IMGUIDED_FILTER: The best delta_ssim to image %d has the size filter [%d %d]. Before SSIM: %.4f. After SSIM: %.4f. Diff = %.2f%% \n', k, best_i, best_i, SP.IP.ssim_wn(SP.pixels_plot(k)), best_ssim_temp, best_delta_ssim_temp*100/SP.IP.ssim_wn(SP.pixels_plot(k)));
                best_ssim_temp = -1; 
                best_delta_ssim_temp = 0;
            end
        
            fclose(fileID);   
end

