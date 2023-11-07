function [roi1, roi2] = get_ref_spectra(SP)

    % imagesc(squeeze(mean((SP.data_processed(1).data_T),[1])));
    % colormap('hot')
    % mask = roipoly;
    % SP = SP.make_raman_spectrum_with_mask(mask);
    % roi1= SP.ramanSpectrumWithMask;
    % 
    % mask2 = roipoly;
    % SP = SP.make_raman_spectrum_with_mask(mask2);
    % roi2= SP.ramanSpectrumWithMask;



    mask_to_int = zeros(50,50);
    spectra_of_pixel = zeros(50,50, 200);

    max_signal_of_image = 0;

    for row = 1:size(mask_to_int, 1)
        for col = 1:size(mask_to_int, 2)
            mask_to_int(row, col) = 1;


            SP = SP.make_raman_spectrum_with_mask(mask_to_int);
            spectra_of_pixel(row, col,:) = SP.ramanSpectrumWithMask(1:200);
            
            if max(spectra_of_pixel(row,col,:)) > max_signal_of_image
                max_signal_of_image = max(spectra_of_pixel(row,col,:));
            end
            mask_to_int(row, col) = 0;

        end
    end
    %Normalize with stronger signal of pixels 
    figure,
    for row = 1:size(mask_to_int, 1)
        for col = 1:size(mask_to_int, 2)
            spectra_of_pixel(row, col,:) = spectra_of_pixel(row, col,:)./max_signal_of_image;
            plot(squeeze(spectra_of_pixel(row, col, :))), hold on
            [row col]
        end
    end

end
