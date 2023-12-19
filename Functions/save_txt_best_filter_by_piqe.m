function save_txt_best_filter_by_piqe(SP)
            cell_with_imgs_wn = SP.IP.mat_img_wn(SP.pixels_plot([1 2 3]));

            %TXT WITH RESULTS
            fileID = fopen(append('results_best_filter_for_piqe_exp',string(SP.xp_number), '.txt'),'w');
            fprintf(fileID, 'Analysis date: %s \n', datetime("now"));
            fprintf(fileID, 'Intern: Paulo Henrique DINIZ FERNANDES\n');
            fprintf(fileID, 'Folder path : %s \n', string(pwd));
            fprintf(fileID, 'Experiment number: %d \n', SP.xp_number);
            fprintf(fileID, 'Tukey Window parameter: %f \n', SP.tukey_window_param);
            fprintf(fileID, 'Percent FWHM: %.4f \n', SP.percent_FWHM);
            fprintf(fileID, 'Window type: %s \n', SP.window2_name);
            fprintf(fileID, 'Ratio Window: %f \n', SP.ratio_window);
            fprintf(fileID, 'Peaks in order: %d, %d, %d, %d, %d \n', SP.pixels_plot(1), SP.pixels_plot(2), SP.pixels_plot(3));
            fprintf(fileID, 'SSIM of peaks without filter: %s, %s, %s \n', num2str(SP.IP.ssim_wn(SP.pixels_plot(1))), num2str(SP.IP.ssim_wn(SP.pixels_plot(2))), num2str(SP.IP.ssim_wn(SP.pixels_plot(3))));
            fprintf(fileID, 'BRISQUE of peaks without filter: %s, %s, %s \n', num2str(SP.IP.brisque_score_wn(SP.pixels_plot(1))), num2str(SP.IP.brisque_score_wn(SP.pixels_plot(2))), num2str(SP.IP.brisque_score_wn(SP.pixels_plot(3))));
            fprintf(fileID, 'NIQE of peaks without filter: %s, %s, %s \n', num2str(SP.IP.niqe_score_wn(SP.pixels_plot(1))), num2str(SP.IP.niqe_score_wn(SP.pixels_plot(2))), num2str(SP.IP.niqe_score_wn(SP.pixels_plot(3))));
            fprintf(fileID, 'PIQE of peaks without filter: %s, %s, %s \n', num2str(SP.IP.piqe_score_wn(SP.pixels_plot(1))), num2str(SP.IP.piqe_score_wn(SP.pixels_plot(2))), num2str(SP.IP.piqe_score_wn(SP.pixels_plot(3))));

            
            fprintf(fileID, '-------------------------------------------\n');

            %GAUSSIEN1
            best_delta_piqe = -1;
            best_i = 1;
            best_sigma = 0;
            for k=1:numel(cell_with_imgs_wn)
                for i=1:10
                    for sigma=0.1:0.05:5
                        img_filtered_temp = SP.IP.gaussian1_filter(cell_with_imgs_wn{k}, [i i], sigma);
                        piqe_temp = piqe(mat2gray(img_filtered_temp));
                        if (SP.IP.piqe_score_wn(SP.pixels_plot(k)) > piqe_temp)
                            delta_piqe_temp = (SP.IP.piqe_score_wn(SP.pixels_plot(k))- piqe_temp);
                            if (delta_piqe_temp > best_delta_piqe)
                                best_delta_piqe = delta_piqe_temp;
                                best_i = i;
                                best_sigma = sigma;
                            end
                        end

                    end
                end
                fprintf(fileID, 'GAUSSIEN1_FILTER: The best delta_piqe to image %d has the size filter [%d %d] and sigma %.2f. Before piqe: %.4f. After piqe: %.4f. Diff = %.2f%% \n', k, best_i, best_i, best_sigma, SP.IP.piqe_score_wn(SP.pixels_plot(k)), best_delta_piqe, (best_delta_piqe-SP.IP.piqe_score_wn(SP.pixels_plot(k)))*100/SP.IP.piqe_score_wn(SP.pixels_plot(k)));
                delta_piqe_temp = -1; 
                best_delta_piqe = 0;
            end
            fclose(fileID);   
end

