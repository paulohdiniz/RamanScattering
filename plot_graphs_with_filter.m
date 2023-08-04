function plot_graphs_with_filter(SP)

    for k = 1:3
        matrix_img_ref = SP.IP.mat_ref;
        temp = SP.IP.mat_img_wn{SP.pixels_plot(k)};
        best_ssim_temp = -1;
        best_i = 1;
        best_j = 1;
            for i=0.1:0.05:3
                for j = 0.1:0.05:3
                    img_filtered_temp = SP.IP.gaussian2_filter(temp, i, j);
                    ssim_temp = ssim(img_filtered_temp,matrix_img_ref);
                    if (ssim_temp > best_ssim_temp)
                        best_ssim_temp = ssim_temp;
                        best_i = i;
                        best_j = j;
                    end
                end
            end
        images(k).image_filtre = ...
        SP.IP.gaussian2_filter(temp, best_i, best_j);
        ssimim(k).imag_filtre = ssim(images(k).image_filtre,matrix_img_ref);
    end


    Total_delay=1/SP.Clock_Freq*SP.N_t*SP.DazzlerTimeConversion;
    time_axis=0:Total_delay/(SP.N_t-1):Total_delay;
    name_of_figure = 'Filtered images';
    h1 = figure('Position', [50 100 900 700], 'Name', name_of_figure);

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
    imagesc(images(1).image_filtre);
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(1))) 'cm^{-1}'] 
        [ 'SSIM(filtered): ' num2str(ssimim(1).imag_filtre)]
        });
    colormap('hot')
    
    subplot(2,3,5);
    imagesc(images(2).image_filtre);
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(2))) 'cm^{-1}'] 
        [ 'SSIM(filtered): ' num2str(ssimim(2).imag_filtre)]
        });
    
    subplot(2,3,6);
    imagesc(images(3).image_filtre);
    xlabel('pixels');
    ylabel('pixels');
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(3))) 'cm^{-1}'] 
        [ 'SSIM(filtered): ' num2str(ssimim(3).imag_filtre)]
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
        append('Deadtime: ', string(SP.deadtime)), ...
        'Filtre: gaussien2', ...
        },...
        'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','bold', ...
        'fontsize',10, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    han.Position(1) = han.Position(1) - abs(han.Position(1) * 0.8); %horizontal indent
    
end

