function plot_similar_pixels_from_rois(SP)

    imagesc(squeeze(mean((SP.data_processed(1).data_T),[1])));
    colormap('hot')
    mask = roipoly;
    SP = SP.make_raman_spectrum_with_mask(mask);
    roi1= SP.ramanSpectrumWithMask;

    mask2 = roipoly;
    SP = SP.make_raman_spectrum_with_mask(mask2);
    roi2= SP.ramanSpectrumWithMask;

    for ii=1:50
        for jj=1:50
            corr_pixel_to_1 = corr(squeeze(abs(SP.hyperspectralRamanImageComplex(1:200,ii,jj))),roi1(1:200));
            pourc_discr_1(ii,jj) = corr_pixel_to_1; 
    
            corr_pixel_to_2 = corr(squeeze(abs(SP.hyperspectralRamanImageComplex(1:200,ii,jj))),roi2(1:200));
            pourc_discr_2(ii,jj) = corr_pixel_to_2; 
            
        end
    end
    figure,imagesc(pourc_discr_1)
    figure,imagesc(pourc_discr_2)

end
