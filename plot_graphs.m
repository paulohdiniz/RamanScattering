function plot_graphs(SP)

    Total_delay=1/SP.Clock_Freq*SP.N_t*SP.DazzlerTimeConversion;
    time_axis=0:Total_delay/(SP.N_t-1):Total_delay;

    figure('Position', [50 100 900 700]);
    subplot(2,3,1);

    hold on,
    plot(SP.data_stitched.t_stitched(:)*1E12,squeeze(mean(SP.data_stitched.data_R,[1 2])).*SP.window2.'); hold off
    xlabel('delay [ps]','fontsize',14);
    ylabel('Raman oscillations','fontsize',14);
    title('Integrated temporal response','fontsize',14);
    text(-0.1,1.1,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
    set(gca,'ytick',[])      

    subplot(2,3,2);
    plot(SP.wn,SP.ramanSpectrum);
    xlabel('Wavenumbers [cm^{-1}]','fontsize',14);
    ylabel('Raman Spectrum','fontsize',14);
    xlim([0,160])
    title('Integrated Raman Spectrum','fontsize',14)
    text(-0.1,1.1,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
    set(gca,'ytick',[]);

    subplot(2,3,3);
    imagesc(squeeze(mean((SP.data_processed(1).data_T),[1])));
    title('DC PD signal','fontsize',14); 
    xlabel('pixels','interpreter','tex');
    ylabel('pixels');
    axis image
    colorbar;
    title('Integrated Raman Spectrum','fontsize',14)
    text(-0.1, 1.1,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)

    subplot(2,3,4);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(1),:,:))));
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title(sprintf('Image at %.1f cm^{-1}',SP.wn(SP.pixels_plot(1))));
    colormap('hot')
    
    subplot(2,3,5);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(2),:,:))));
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title(sprintf('Image at %.1f cm^{-1}',SP.wn(SP.pixels_plot(2))));
    
    subplot(2,3,6);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(3),:,:))));
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title(sprintf('Image at %.1f cm^{-1}',SP.wn(SP.pixels_plot(3))));


end

