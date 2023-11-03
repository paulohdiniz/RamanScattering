function plot_best_ssim_by_ratio_window(SPtemp)

    name_of_figure = 'Best SSIM x RW';
    h1 = figure('Position', [50 100 1000 1300], 'Name', name_of_figure);

    best_ssim =[]; %yAxis
    step = 0.1; 
    ratios = step:step:1;
    %windows = {'barthannwin','blackman','flattopwin','kaiser','rectwin','ones'};
    %colors = {'red','green','blue',"#7E2F8E","cyan", "#000000"};
    windows = {SPtemp.window2_name};
    for j = 1:length(windows)
        for i = 1:numel(ratios)
            SP(i) = SPtemp.copy(); %copying the data internally is much faster than loading all the files every time
            SP(i).window2_name = windows{j};
            SP(i).ratio_window = ratios(i);
            SP(i) = SP(i).window_overlap_to_test(SP(i).tukey_window_param,SP(i).pourc_pulse_width);
            SP(i) = SP(i).Tnorm_and_center_data(1,0);
            SP(i) = SP(i).stitch_time_axis_T_with_interp(SP(i).interp_method);
            SP(i) = SP(i).pick_fourier_window(SP(i).window2_name); 
            SP(i) = SP(i).FT(SP(i).data_stitched.t_stitched, permute(SP(i).data_stitched.data_R,[3 1 2]).*repmat(SP(i).window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman SP(i)ectrum is arbitrary units
            SP(i) = SP(i).make_raman_spectrum();
            SP(i) = SP(i).calculated_images_scores_per_wn();
            best_ssim(i) = (SP(i).IP.peaks_ssim(1)+SP(i).IP.peaks_ssim(2)+SP(i).IP.peaks_ssim(3))/3;
            
            values_to_plot = 0.1:0.1:1;
            if ismember(ratios(i), values_to_plot)
                [((j==1)*(1.0 + step -ratios(i)) + (j~=1)*0) ((j==2)*(1.0 + step -ratios(i)) + (j~=2)*0) ((j==3)*(1.0 + step -ratios(i)) + (j~=3)*0)]
                subplot(2,1,2);
                plot(SP(i).wn, SP(i).ramanSpectrum, ...
                    'DisplayName', append(SP(i).window2_name,':', string(SP(i).ratio_window)),...
                    'Color',[((j==1)*(1.0 + step -ratios(i)) + (j~=1)*0) ((j==2)*(1.0 + step -ratios(i)) + (j~=2)*0) ((j==3)*(1.0 + step -ratios(i)) + (j~=3)*0)], ...
                    LineWidth=0.5 + j*0.5); hold on
                xlabel('Wavenumbers [cm^{-1}]','fontsize',14);
                ylabel('Raman Spectrum','fontsize',14);
                xlim([0,160])
                legend();
            end


        end
        subplot(2,1,1);
        plot(ratios, best_ssim, 'DisplayName', SP(i).window2_name); hold on
        
    end
    xlabel('Ratio Window','fontsize',14);
    ylabel('Avg top 3 SSIMs','fontsize',14);
    legend();

        %Putting Parameters
    han=axes(h1,'visible','off');
    han.YLabel.Visible='on';
    ylabel(han,{ ...
        append('Exp: ', string(SPtemp.xp_number)), ...
        SPtemp.function_generator, ...
        SPtemp.lockin_parameters, ...
        append('Tukey ratio: ', string(SPtemp.tukey_window_param)), ...
        append('Pourc pulse: ', string(SP.pourc_pulse_width))}, ...
        'Rotation',0, ...
        'interpreter','none', ...
        'fontweight','bold', ...
        'fontsize',9, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    han.Position(1) = han.Position(1) - abs(han.Position(1) * 0.8); %horizontal indent

end

