function plot_graphs_with_ifft(SP)
    
    if isprop(SP,'xp_number')
        name_of_figure = append('Exp : ', string(SP.xp_number));
    else
        name_of_figure = "Plot Graphs with ifft";
    end

    h1 = figure('Position', [50 100 1300 800],'Name', name_of_figure);
    
    subplot(6,4,[1,2,5,6,9,10]);

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

    subplot(6,4,[13, 14, 17, 18, 21, 22]);
    plot(SP.wn,SP.ramanSpectrum);
    hold on,plot(SP.wns_plot,SP.ramanSpectrum(SP.pixels_plot),'ro')
    xlabel('Wavenumbers [cm^{-1}]','fontsize',8);
    ylabel('Raman Spectrum','fontsize',8);
    xlim([0,160])
    title('Integrated Raman Spectrum','fontsize',8)
    text(-0.1,1.1,'','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',8)
    set(gca,'ytick',[]);

    time_half_axis = SP.data_stitched.t_stitched(1:floor(length(SP.data_stitched.t_stitched)/5)).*1E12;
    
    if isempty(SP.signalIFFT)
        SP = SP.make_signal_ifft();
    end
    subplot(6,4, [3,4,7,8])
    plot(time_half_axis, squeeze(mean(real(SP.signalIFFT(1).signal(1:length(time_half_axis),:,:)),[2 3])))

    subplot(6,4, [11,12,15,16])
    plot(time_half_axis, squeeze(mean(real(SP.signalIFFT(2).signal(1:length(time_half_axis),:,:)),[2 3])))

    subplot(6,4, [19,20,23,24])
    plot(time_half_axis, squeeze(mean(real(SP.signalIFFT(3).signal(1:length(time_half_axis),:,:)),[2 3])))
    xlabel('ps','fontsize',10);

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
        'fontsize',9, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    han.Position(1) = han.Position(1) - abs(han.Position(1) * 0.8); %horizontal indent

end

