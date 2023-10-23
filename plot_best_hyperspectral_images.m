function plot_best_hyperspectral_images(SP)

name_of_figure = append('Exp : ', string(SP.xp_number));

[t_values_for_stop, t_subsections, y_subsections, t, ampli, freq_of_peaks] = get_ideal_times_by_spectrogram(SP);

for i=1:numel(t_values_for_stop)
    SPs(i) = SP.copy(); %copying the data internally is much faster than loading all the files every time
    SPs(i).ratio_window = t_values_for_stop(i)/max(SPs(i).data_stitched.t_stitched);
    SPs(i) = SPs(i).window_overlap_to_test(SPs(i).tukey_window_param,SPs(i).deadtime);
    SPs(i) = SPs(i).Tnorm_and_center_data(1,0,SPs(i).deadtime);
    SPs(i) = SPs(i).stitch_time_axis_T_with_interp(SPs(i).interp_method);
    SPs(i) = SPs(i).pick_fourier_window(SPs(i).window2_name); 
    SPs(i) = SPs(i).FT(SPs(i).data_stitched.t_stitched, permute(SPs(i).data_stitched.data_R,[3 1 2]).*repmat(SPs(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
    SPs(i) = SPs(i).make_raman_spectrum();
    SPs(i) = SPs(i).points_to_plot_by_frequency();   
end
h1 = figure('Position', [50 100 1300 800],'Name', name_of_figure);
for i=1:numel(t_values_for_stop)

    signal = squeeze(mean(SPs(i).data_stitched.data_R,[1 2]));
    signal_norm=signal./max(signal(:));
    signal_windowed = signal_norm.*SPs(i).window2.';

    %Pos 4*i -3
    subplot(3,4,4*i -3);
    hold on, plot(SPs(i).data_stitched.t_stitched(:)*1E12,signal_windowed, 'Color', 'blue',LineWidth=1.5);
    hold on, plot(SPs(i).data_stitched.t_stitched(:)*1E12,signal_norm,"--", 'Color', '#006400',LineWidth=1);
    hold on, plot(SPs(i).data_stitched.t_stitched(:)*1E12, SPs(i).window2, "--",'Color', 'red',LineWidth=2);
    hold off
    xlh = xlabel({append('Time opt: ', sprintf('%.2fps', t_values_for_stop(i)*1E12)),'delay [ps]'},'fontsize',8);
    %xlh.Position(1) = xlh.Position(1) + abs(xlh.Position(1) * 0.95); %TODO:regler les distances
    ylabel('Raman oscillations','fontsize',8);
    title('Integrated temporal response','fontsize',8);
    text(-0.1,1.1,'','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',8)
    set(gca,'ytick',[])

    %Pos 4*i -2
    subplot(3,4,4*i -2);
    plot(SPs(i).wn,SPs(i).ramanSpectrum);
    hold on,plot(SP.wns_plot(i),SPs(i).ramanSpectrum(SP.pixels_plot(i)),'ro')
    xlabel('Wavenumbers [cm^{-1}]','fontsize',8);
    ylabel('Raman Spectrum','fontsize',8);
    xlim([0,160])
    title('Integrated Raman Spectrum','fontsize',8)
    text(-0.1,1.1,'','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',8)
    set(gca,'ytick',[]);

    %Pos 4*i -1
    subplot(3,4,4*i -1);
    plot(t*1E12, ampli{i}), hold on,
    ylabel('Amplitude [u.a]','fontsize',8)
    xlabel({append('Freq: ', sprintf('%.2f cm^{-1}', freq_of_peaks(i)/100/SP.c)),'delay [ps]'},'fontsize',8);
    plot(t_subsections{i}*1E12, y_subsections{i}), hold on,
    area(t_subsections{i}*1E12, y_subsections{i}, 'FaceColor', 'b'), hold on
    % if(t_values_for_stop(i) )
    % end
    area(t(find(t >= t_values_for_stop(i)))*1E12, ampli{i}(find(t >= t_values_for_stop(i)))', 'FaceColor', 'r')
    
    %Pos 4*i
    subplot(3,4,4*i)
    imagesc(squeeze(abs(SPs(i).hyperspectralRamanImageComplex(SP.pixels_plot(i),:,:))));
    axis image
    xlabel('pixels','fontsize',6);
    ylabel('pixels','fontsize',6);
    colorbar;
    title({
        ['Image at ' num2str(SPs(i).wn(SP.pixels_plot(i))) 'cm^{-1}'] 
        [ 'SSIM: ' num2str(SPs(i).IP.ssim_wn(SP.pixels_plot(i)))]
        },'fontsize',8);
    colormap('hot')

end

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
    'fontsize',9, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom');
han.Position(1) = han.Position(1) - abs(han.Position(1) * 0.8); %horizontal indent
end

