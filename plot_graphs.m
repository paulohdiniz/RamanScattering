function plot_graphs(SP)

    name_of_figure = append('Exp : ', string(SP.xp_number));
    h1 = figure('Position', [50 100 1300 800],'Name', name_of_figure);
    
    subplot(6,3,[1,2,4,5,7,8]);

    signal = squeeze(mean(SP.data_stitched.data_R,[1 2]));
    signal_norm=signal./max(signal(:));

    signal_windowed = signal_norm.*SP.window2.';
    
    hold on, plot(SP.data_stitched.t_stitched(:)*1E12,signal_windowed, 'Color', 'blue',LineWidth=1.5);
    hold on, plot(SP.data_stitched.t_stitched(:)*1E12,signal_norm,"--", 'Color', '#006400',LineWidth=1);
    hold on, plot(SP.data_stitched.t_stitched(:)*1E12, SP.window2, "--",'Color', 'red',LineWidth=2);
    hold off

    xlh = xlabel('delay [ps]','fontsize',8);
    xlh.Position(1) = xlh.Position(1) + abs(xlh.Position(1) * 0.95);
    ylabel('Raman oscillations','fontsize',8);
    title('Integrated temporal response','fontsize',8);
    text(-0.1,1.1,'','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',8)
    set(gca,'ytick',[])      

    subplot(6,3,[10,11,13,14,16,17]);
    plot(SP.wn,SP.ramanSpectrum);
    hold on,plot(SP.wns_plot,SP.ramanSpectrum(SP.pixels_plot),'ro')
    xlabel('Wavenumbers [cm^{-1}]','fontsize',8);
    ylabel('Raman Spectrum','fontsize',8);
    xlim([0,160])
    title('Integrated Raman Spectrum','fontsize',8)
    text(-0.1,1.1,'','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',8)
    set(gca,'ytick',[]);

    subplot(6,3,[3,6]);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(1),:,:))));
    axis image
    xlabel('pixels','fontsize',6);
    ylabel('pixels','fontsize',6);
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(1))) 'cm^{-1}'] 
        [ 'SSIM: ' num2str(SP.IP.ssim_wn(SP.pixels_plot(1)))]
        },'fontsize',8);
    colormap('hot')
    
    subplot(6,3,[9,12]);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(2),:,:))));
    axis image
    xlabel('pixels','fontsize',6);
    ylabel('pixels','fontsize',6);
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(2))) 'cm^{-1}'] 
        [ 'SSIM: ' num2str(SP.IP.ssim_wn(SP.pixels_plot(2)))]
        },'fontsize',8);
    
    subplot(6,3,[15,18]);
    imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(3),:,:))));
    axis image
    xlabel('pixels','fontsize',6);
    ylabel('pixels','fontsize',6);
    colorbar;
    title({
        ['Image at ' num2str(SP.wn(SP.pixels_plot(3))) 'cm^{-1}'] 
        [ 'SSIM: ' num2str(SP.IP.ssim_wn(SP.pixels_plot(3)))]
        },'fontsize',8);
    
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
        append('ScoreCriteria1: ', string(SP.scoreCriteria(1))), ...
        append('ScoreCriteria2: ', string(SP.scoreCriteria(2))), ...
        append('ScoreCriteria3: ', string(SP.scoreCriteria(3))), ...
        append('ScoreCriteria4: ', string(SP.scoreCriteria(4))), ...
        append('ScoreCriteria5: ', string(SP.scoreCriteria(5))) ...
        }, ...
        'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','bold', ...
        'fontsize',9, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    han.Position(1) = han.Position(1) - abs(han.Position(1) * 0.8); %horizontal indent
end

