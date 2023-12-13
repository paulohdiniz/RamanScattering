function plot_spectrogram(SP)

    %Default Values of Spectrogram
    temp_t=SP.data_stitched.t_stitched;
    temp_data=squeeze(mean(SP.data_stitched.data_R,[1 2]));
    freqLim = 150; %in cm^-1

    M = round(numel(temp_t)*0.5); %not optimized yet
    L = M - 5 ;    %not optimized yet, PS: Changing L only changes the qtd of points of t, s and f not changing its size
    Ndft = max(512,2^nextpow2(M));
    g = blackman(M); 
    [s,f,t, p] = spectrogram(temp_data,g,L,Ndft,SP.Fs);
    [r,~,~,q] = spectrogram(temp_data,g,L,Ndft,SP.Fs,"power");

    f_in_cm =f/100/SP.c;
    t_in_ps = t.*1E12;

    name_of_figure = append('Exp : ', string(SP.xp_number));
    h1 = figure('Position', [50 100 1400 900],'Name', name_of_figure);

    subplot(4,9,[1,2,3])
        pspectrum(temp_data,temp_t, Leakage=0.5,FrequencyLimits=[0,freqLim*100*SP.c])
        subtitle("Leakage: 0.5")
    subplot(4,9,[4,5,6])
        pspectrum(temp_data,temp_t, Leakage=0.75,FrequencyLimits=[0,freqLim*100*SP.c])
        subtitle("Leakage: 0.75")
        ylabel('')
    subplot(4,9,[7,8,9])
        pspectrum(temp_data,temp_t, Leakage=1,FrequencyLimits=[0,freqLim*100*SP.c])
        ylabel('')
        subtitle("Leakage: 1")
    %Pos1
    subplot(4,9,[10,11,12,13])
        pspectrum(temp_data,SP.Fs,"spectrogram", ...
        TimeResolution=M/SP.Fs,OverlapPercent=L/M*100, ...
        Leakage=0.7 ...
        ,FrequencyLimits=[0,freqLim*100*SP.c])
        title("pspectrum")
        %ylabel('Wavenumbers [cm^{-1}]','fontsize',8);
        axis image

    subplot(4,9,[15,16,17,18])
        pspectrum(temp_data,SP.Fs,"persistence", ...
        TimeResolution=M/SP.Fs,OverlapPercent=L/M*100, ...
        Leakage=0.7 ...
        ,FrequencyLimits=[0,freqLim*100*SP.c])
        title("Persistence spectrum")
        %xlabel('Wavenumbers [cm^{-1}]','fontsize',8);
    %Pos2
    subplot(4,9,[19,20,21,22,28,29,30,31])
        waterplot(s,f_in_cm,t_in_ps)
    subplot(4,9,[24,25,26,27,33,34,35,36])
        %Compute Segment PSDs and Power Spectra
        waterfall(f_in_cm,t_in_ps,abs(s)'.^2)
        set(gca,XScale="log",...
        XDir="reverse",View=[30 50], XLim=[0 freqLim] )
        title('Power spectrum','fontsize',8)
        xlabel('Wavenumbers [cm^{-1}]','fontsize',8);
        ylabel('delay [ps]','fontsize',8)
        zlabel('????? [u.a]','fontsize',8)

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

