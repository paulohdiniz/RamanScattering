function plot_spectrogram(SP)

    %Default Values of Spectrogram
    temp_t=SP.data_stitched.t_stitched;
    temp_data=squeeze(mean(SP.data_stitched.data_R,[1 2]));

    freqLim = 150; %in cm^-1

    M = round(numel(temp_t)*0.5); %pas optimise encore
    L = M - 5 ;    %pas optimise encore, PS: Changer L change seulement la qtd de points de t, s et f ne chageant pas sa taille
    Ndft = max(512,2^nextpow2(M));
    g = blackman(M); 
    [s,f,t, p] = spectrogram(temp_data,g,L,Ndft,SP.Fs);
    [r,~,~,q] = spectrogram(temp_data,g,L,Ndft,SP.Fs,"power");

    f_in_cm =f/100/SP.c;
    t_in_ps = t.*1E12;

    name_of_figure = append('Exp : ', string(SP.xp_number));
    h1 = figure('Position', [50 100 500 800],'Name', name_of_figure);
    %Pos1
    subplot(3,4,[1,2,3,4])
        pspectrum(temp_data,SP.Fs,"spectrogram", ...
        TimeResolution=M/SP.Fs,OverlapPercent=L/M*100, ...
        Leakage=0.7 ...
        ,FrequencyLimits=[0,freqLim*100*SP.c])
        title("pspectrum")
        axis image
    %Pos2
    subplot(3,4,[5,6,7,8])
        waterplot(s,f_in_cm,t_in_ps)
    subplot(3,4,[9,10,11,12])
        %Compute Segment PSDs and Power Spectra
        waterfall(f_in_cm,t_in_ps,abs(r)'.^2)
        set(gca,XScale="log",...
        XDir="reverse",View=[30 50], XLim=[0 freqLim] )
        title('Power spectrum???','fontsize',8)
        xlabel('Wavenumbers [cm^{-1}]','fontsize',8);
        ylabel('delay [ps]','fontsize',8)
        zlabel('????? [u.a]','fontsize',8)

    % 
    % for M=31:10:800
    %     for L=M-30:M-10
    %         if (L < M)
    %             g = blackman(M); 
    % 
    %             Ndft = max(256,2^nextpow2(M));
    %             texto = [M L Ndft];
    %             [s,f,t] = spectrogram(temp_data,g,L,Ndft,SP.Fs);
    %             waterplot(s,f,t,texto),            
    %             pause(0.001)
    %             % close gcf,
    %         end
    %     end
    % end

    % Waterfall plot of spectrogram
    function waterplot(s,f,t)
        x =waterfall(f,t,abs(s)'.^2);
        set(gca,XDir="reverse",View=[30 50], XLim=[0 freqLim] )
        title(append('Spectrogram. ','M=',string(M),', ', 'L=', string(L), ', ', 'Ndft=', string(Ndft) ),'fontsize',8)
        xlabel('Wavenumbers [cm^{-1}]','fontsize',8);
        ylabel('delay [ps]','fontsize',8)
        zlabel('Amplitude [u.a]','fontsize',8)
    end
end

