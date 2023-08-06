function save_graphs_as_PDF(SP)

    Total_delay=1/SP.Clock_Freq*SP.N_t*SP.DazzlerTimeConversion;
    time_axis=0:Total_delay/(SP.N_t-1):Total_delay;
    name_of_figure = append('Exp : ', string(SP.xp_number), '\n');
    h1 = figure('Position', [50 100 900 700],'visible','off', 'Name', name_of_figure);
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
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(1))) 'cm^{-1}'] 
        [ 'SSIM: ' num2str(SP.IP.peaks_ssim(1))]
        });
    colormap('hot')
    
    subplot(2,3,5);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(2),:,:))));
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(2))) 'cm^{-1}'] 
        [ 'SSIM: ' num2str(SP.IP.peaks_ssim(2))]
        });
    
    subplot(2,3,6);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(3),:,:))));
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(3))) 'cm^{-1}'] 
        [ 'SSIM: ' num2str(SP.IP.peaks_ssim(3))]
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
        append('Deadtime: ', string(SP.deadtime))}, ...
        'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','bold', ...
        'fontsize',8, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    han.Position(1) = han.Position(1) - abs(han.Position(1) * 0.8); %horizontal indent
    
    if exist('test.pdf', 'file')
            exportgraphics(h1,'test.pdf',"Append",true, 'BackgroundColor','none');
    else
            exportgraphics(h1,'test.pdf', 'BackgroundColor','none');
            %print(h1,'test.pdf','-dpdf', '-bestfit')
    end

end
