classdef Image_Processing
    
    properties
        img_ref
        mat_ref

        img_wn
        mat_img_wn

        ssim_wn
        mean_ssim_wn
        peaks_ssim_wn
        peaks_ssim

        imag_filtered

    end
    
    methods
        %AVG FILTER
        function [imag_filtered] = average_filter(~, imag_raw, hsize) % >1 ok!
            imag_filtered = imfilter(imag_raw,fspecial('average', hsize));
        end
        function [best_img_filtered] = get_best_img_avg_filter(~,SP, imag_raw)           
            best_ssim_temp = -1; 
            best_i = 1;
            best_j = 1;
            for i=1:25
                for j = 1:25
                    img_filtered_temp = SP.IP.average_filter(imag_raw, [i j]);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_i = i;
                        best_j = j;
                        best_img_filtered = img_filtered_temp;
                    end
                end
            end
        end


        function [imag_filtered] = gaussian1_filter(~, imag_raw, hsize,sigma) % hsize > 1 et sigma > 0 ok!
            imag_filtered = imfilter(imag_raw,fspecial('gaussian', hsize, sigma));
        end
        function [imag_filtered] = disk_filter(~, imag_raw, radius) % radius > 1 et < 12 ok!
            imag_filtered = imfilter(imag_raw,fspecial('disk', radius));
        end
        function [imag_filtered] = laplacian_filter(~, imag_raw, alpha) % alpha [0,1] ok!
            imag_filtered = imfilter(imag_raw,fspecial('laplacian', alpha));
        end
        function [imag_filtered] = log_filter(~, imag_raw, hsize, sigma) %hsize > 1 et sigma > 0 ok!
            imag_filtered = imfilter(imag_raw,fspecial('log', hsize, sigma));
        end
        function [imag_filtered] = motion_filter(~, imag_raw, len, theta) % len > 1 , 0 < theta < 360 ok!
            imag_filtered = imfilter(imag_raw,fspecial('motion', len, theta));
        end
        function [imag_filtered] = prewitt_filter(~, imag_raw) %without parameters ok
            imag_filtered = imfilter(imag_raw,fspecial('prewitt'));
        end
        function [imag_filtered] = sobel_filter(~, imag_raw) %without parameters ok
            imag_filtered = imfilter(imag_raw,fspecial('sobel'));
        end

        function [imag_filtered] = medfilt2_filter(~, imag_raw , m, n) %ok
            imag_filtered = medfilt2(imag_raw,[m n]);
        end
        function [imag_filtered] = wiener2_filter(~, imag_raw, m, n) %ok
            [imag_filtered,~] = wiener2(imag_raw,[m n]);
        end

        function [imag_filtered] = gaussian2_filter(~, imag_raw, sigma_m ,sigma_n ) %sigma > 0 OK!
            imag_filtered = imgaussfilt(imag_raw, [sigma_m sigma_n]);
        end
        function [best_img_filtered] = get_best_img_gauss2_filter(~,SP, imag_raw)  
            best_ssim_temp = -1;
            for i=0.1:0.05:3
                for j = 0.1:0.05:3
                    img_filtered_temp = SP.IP.gaussian2_filter(imag_raw, i, j);
                    ssim_temp = ssim(img_filtered_temp,SP.IP.mat_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_img_filtered = img_filtered_temp;
                    end
                end
            end
        end

        function [imag_filtered] = mode_filter(~, imag_raw, size) %size = vector of positive odd integers
            imag_filtered = modefilt(imag_raw, size); %OK !
        end
        function [imag_filtered] = ordfilt2_filter(~, imag_raw, order) % 2 <= order <= 9 ok!
            imag_filtered = ordfilt2(imag_raw, order, ones(3,3));
        end
        function [imag_filtered] = std_filter(~, imag_raw, n) %n neighborhood ok!
            imag_filtered = stdfilt(imag_raw, n);
        end
        function [imag_filtered] = entropy_filter(~, imag_raw, n) %n neighborhood ok!
            imag_filtered = entropyfilt(imag_raw, n);
        end
        function [imag_filtered] = range_filter(~, imag_raw, n) %n neighborhood ok!
            imag_filtered = rangefilt(imag_raw, n);
        end

        function [imag_filtered] = imguided_filter(~, imag_raw, m , n) % 2 <= N
            imag_filtered = imguidedfilter(imag_raw, NeighborhoodSize=[m n]);
        end

        % The image size should be greater than or equal to 64-by-64.
        function [imag_filtered] = imdiffuse_filter(~, imag_raw) %without parameters
            [gradThresh,numIter] = imdiffuseest(imag_raw,'ConductionMethod','quadratic');
            imag_filtered = imdiffusefilt(imag_raw,'ConductionMethod','quadratic', ...
            'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
        end

    end

end

