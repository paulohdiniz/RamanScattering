function plot_graphs_with_roi(SP)

    Total_delay=1/SP.Clock_Freq*SP.N_t*SP.DazzlerTimeConversion;
    time_axis=0:Total_delay/(SP.N_t-1):Total_delay;

    h1 = figure('Position', [25 100 975 800]);
    subplot(2,3,1);

    hold on,
    plot(SP.data_stitched.t_stitched(:)*1E12,squeeze(mean(SP.data_stitched.data_R,[1 2])).*SP.window2.'); hold off
    xlabel('delay [ps]','fontsize',14);
    ylabel('Raman oscillations','fontsize',14);
    title('Integrated temporal response','fontsize',14);
    text(-0.1,1.1,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
    set(gca,'ytick',[])      

    subplot(2,3,2);
    plot(SP.wn,SP.ramanSpectrum, 'DisplayName', 'Total Spectrum', 'Color', 'blue');
    hold on,plot(SP.wns_plot,SP.ramanSpectrum(SP.pixels_plot),'ro', 'HandleVisibility','off')
    legend();
    xlabel('Wavenumbers [cm^{-1}]','fontsize',14);
    ylabel('Raman Spectrum','fontsize',14);
    xlim([0,160])
    title('Integrated Raman Spectrum','fontsize',14)
    text(-0.1,1.1,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
    axisRamSpec = gca;
    %set(gca,'ytick',[]);

    subplot(2,3,3);
    imagesc(squeeze(mean((SP.data_processed(1).data_T),[1])));
    title('DC PD signal','fontsize',14); 
    axis image
    xlabel('pixels','interpreter','tex');
    ylabel('pixels');
    axis image
    colorbar;
    title('Integrated Raman Spectrum','fontsize',14)
    text(-0.1, 1.1,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
    figure_DATAT = gcf;

    subplot(2,3,4);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(1),:,:))));
    axis image
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(1))) 'cm^{-1}'] 
        [ 'SSIM: ' num2str(SP.IP.ssim_wn(SP.pixels_plot(1)))]
        });
    colormap('hot')
    
    subplot(2,3,5);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(2),:,:))));
    axis image
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(2))) 'cm^{-1}'] 
        [ 'SSIM: ' num2str(SP.IP.ssim_wn(SP.pixels_plot(2)))]
        });
    
    subplot(2,3,6);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(3),:,:))));
    axis image
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(3))) 'cm^{-1}'] 
        [ 'SSIM: ' num2str(SP.IP.ssim_wn(SP.pixels_plot(3)))]
        });

    %Putting Parameters
    han=axes(h1,'visible','off'); 
    han.YLabel.Visible='on';
    ylabel(han,{ ...
        append('Exp: ', string(SP.xp_number)), ...
        SP.function_generator, ...
        SP.lockin_parameters, ...
        append('Window: ', string(SP.window2_name)) ...
        append('Ratio window: ', string(SP.ratio_window)), ...
        append('Tukey ratio: ', string(SP.tukey_window_param)), ...
        append('Perc FWHM: ', string(SP.percent_FWHM))}, ...
        'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','bold', ...
        'fontsize',10, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    han.Position(1) = han.Position(1) - abs(han.Position(1) * 0.83); %horizontal indent
    han.Position(2) = han.Position(2) - abs(han.Position(2) * 1); %vertical indent

    figure('Name','Select the desired pixel polygon:', 'Position', [1000 500 400 400]);
    imagesc(squeeze(mean((SP.data_processed(1).data_T),[1])));
    colormap('hot')

    % ROI 1
    mask = roipoly;
    SP = SP.make_raman_spectrum_with_mask(mask);
    tempmask1= SP.ramanSpectrumWithMask;
    plot(axisRamSpec, SP.wn,tempmask1, "-", 'DisplayName','ROI Spectrum','Color', '#006400');

    % ROI 2
    mask2 = roipoly;
    SP = SP.make_raman_spectrum_with_mask(mask2);
    tempmask2= SP.ramanSpectrumWithMask;
    plot(axisRamSpec, SP.wn,tempmask2, "-",'DisplayName','ROI 2 Spectrum','Color', 'red');

    figure ('Name', 'Raman Spectrum with mask selected','Position',[1000 100 700 500]),
    plot(SP.wn,SP.ramanSpectrum, 'DisplayName','Total Spectrum', 'Color', 'blue', LineWidth=2);
    hold on,
    plot(SP.wn,tempmask1, "-",'DisplayName','ROI Spectrum', 'Color', '#006400',LineWidth=2);
    hold on,
    plot(SP.wn,tempmask2, "-", 'DisplayName','ROI 2 Spectrum','Color', 'red',LineWidth=2);
    legend();
    xlabel('Wavenumbers [cm^{-1}]','fontsize',14);
    ylabel('Raman Spectrum','fontsize',14);
    xlim([0,160])
    title('Integrated Raman Spectrum with mask selected','fontsize',14)
    text(-0.1,1.1,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
    set(gca,'ytick',[]);
    
end