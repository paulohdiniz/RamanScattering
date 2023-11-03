function plot_graphs_with_transmission(SP)

    name_of_figure = ' ';
    h2 = figure('Position', [50 100 1300 1000], 'Name', name_of_figure);

    subplot(2,3,1);

    signal = squeeze(mean(SP.data_stitched.data_R,[1 2]));
    signal_norm=signal./max(signal(:));

    signal_windowed = signal_norm.*SP.window2.';
    
    hold on, plot(SP.data_stitched.t_stitched(:)*1E12,signal_windowed, 'Color', 'blue',LineWidth=1);
    hold on, plot(SP.data_stitched.t_stitched(:)*1E12,signal_norm,"--", 'Color', '#006400',LineWidth=1);
    hold on, plot(SP.data_stitched.t_stitched(:)*1E12, SP.window2, "--",'Color', 'red',LineWidth=2);
    hold off

    xlabel('delay [ps]','fontsize',14);
    ylabel('Raman oscillations','fontsize',14);
    title('Integrated temporal response','fontsize',14);
    text(-0.1,1.1,'','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
    set(gca,'ytick',[])      

    subplot(2,3,2);
    plot(SP.wn,SP.ramanSpectrum);
    hold on,plot(SP.wns_plot,SP.ramanSpectrum(SP.pixels_plot),'ro')
    %hold on, plot(SP.wns_plot, SP.FT())
    xlabel('Wavenumbers [cm^{-1}]','fontsize',14);
    ylabel('Raman Spectrum','fontsize',14);
    xlim([0,160])
    title('Integrated Raman Spectrum','fontsize',14)
    text(-0.1,1.1,'','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
    %set(gca,'ytick',[]);

    subplot(2,3,3);
    imagesc(squeeze(mean((SP.data_processed(1).data_T),[1])));
    axis image
    title('DC PD signal','fontsize',14); 
    xlabel('pixels','interpreter','tex');
    ylabel('pixels');
    colorbar;
    title('Transmission Image','fontsize',14)
    text(-0.1, 1.1,'','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)

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
    han=axes(h2,'visible','off'); 
    han.YLabel.Visible='on';
    ylabel(han,{ ...
        append('Exp: ', string(SP.xp_number)), ...
        SP.function_generator, ...
        SP.lockin_parameters, ...
        append('Window: ', string(SP.window2_name)) ...
        append('Ratio window: ', string(SP.ratio_window)), ...
        append('Tukey ratio: ', string(SP.tukey_window_param)), ...
        append('Pourc pulse: ', string(SP.pourc_pulse_width)), ...
        },...
        'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','bold', ...
        'fontsize',10, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    han.Position(1) = han.Position(1) - abs(han.Position(1) * 0.85); %horizontal indent
    han.Position(2) = han.Position(2) - abs(han.Position(2) * 1); %vertical indent
end

