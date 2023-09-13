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

        %No-references image quality score
        piqe_score_wn %https://fr.mathworks.com/help/images/ref/piqe.html
        brisque_score_wn %https://fr.mathworks.com/help/images/ref/brisque.html
        niqe_score_wn %https://fr.mathworks.com/help/images/ref/niqe.html


    end
    
    methods
        %AVG FILTER
        function [imag_filtered] = average_filter(~, imag_raw, hsize) % >1 
            imag_filtered = imfilter(imag_raw,fspecial('average', hsize));
        end
        function [best_img_filtered] = get_best_img_avg_filter(~,SP, imag_raw)           
            best_ssim_temp = -1; 
            for i = 1:10
                img_filtered_temp = SP.IP.average_filter(imag_raw, [i i]);
                ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                if (ssim_temp > best_ssim_temp)
                    best_ssim_temp = ssim_temp;
                    best_img_filtered = img_filtered_temp;
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

        function [imag_filtered] = medfilt2_filter(~, imag_raw , m, n) 
            imag_filtered = medfilt2(imag_raw,[m n]);
        end

        function [imag_filtered] = wiener2_filter(~, imag_raw, m, n) 
            [imag_filtered,~] = wiener2(imag_raw,[m n]);
        end

        function [imag_filtered] = gaussian2_filter(~, imag_raw, sigma_m ,sigma_n ) %sigma > 0 OK!
            imag_filtered = imgaussfilt(imag_raw, [sigma_m sigma_n]);
        end

        function [best_img_filtered] = get_best_img_gauss2_filter(~,SP, imag_raw)  
            best_ssim_temp = -1;
            for i = 0.1:0.05:5
                img_filtered_temp = SP.IP.gaussian2_filter(imag_raw, i, i);
                ssim_temp = ssim(img_filtered_temp,SP.IP.mat_ref);
                if (ssim_temp > best_ssim_temp)
                    best_ssim_temp = ssim_temp;
                    best_img_filtered = img_filtered_temp;
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

        function [imag_filtered] = imguided_filter(~, imag_raw, m , n) % 2 <= N
            imag_filtered = imguidedfilter(imag_raw, NeighborhoodSize=[m n]);
        end

        % The image size should be greater than or equal to 64-by-64.
        function [imag_filtered] = imdiffuse_filter(~, imag_raw) %without parameters
            [gradThresh,numIter] = imdiffuseest(imag_raw,'ConductionMethod','quadratic');
            imag_filtered = imdiffusefilt(imag_raw,'ConductionMethod','quadratic', ...
            'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
        end

        function [rgb] = imagesc2rgb(~,im)

            h = imagesc(im); % imagesc handle
            cdata = h.CData; % get image data (if you don't have variable 'im')
            cm = colormap(h.Parent); % get axes colormap
            n = size(cm,1); % number of colors in colormap
            c = linspace(h.Parent.CLim(1),h.Parent.CLim(2),n); % intensity range
            ind = reshape(interp1(c,1:n,im(:),'nearest'),size(im)); % indexed image
            rgb = ind2rgb(ind,cm); % rgb image
        
        end

        function piqe_score = get_piqe_score(~, img)
            piqe_score = piqe(img);
        end

    end

end

