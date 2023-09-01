function plot_spectrogram(SP)

    %Default Values of Spectrogram
    temp_t=SP.data_stitched.t_stitched;
    temp_data=squeeze(mean(SP.data_stitched.data_R,[1 2]));
    %temp_data=temp_data./max(temp_data(:));
    %figure,plot(temp_t,temp_data)
    % 
    % temp=SP.ramanSpectrum;
    % temp=temp./max(temp(:));
    % 
    % wn_temp=SP.wn(SP.wn<SP.Max_wn_to_search_peaks);
    % temp=(temp(SP.wn<SP.Max_wn_to_search_peaks));
    % 
    % wn_temp2=wn_temp(wn_temp>SP.Min_wn_to_search_peaks);
    % temp=(temp(wn_temp>SP.Min_wn_to_search_peaks));
    
    % [s,w,t]=spectrogram(temp_data,'yaxis');
    % figure,imagesc(w,t,10*log10(abs(s)).')
    % figure,waterplot(s,w,t)
    M = 540; 
    L = 520 ;
    Ndft = 1024;
    Ndft = max(512,2^nextpow2(M));
    g = blackman(M); 
    text = [M L Ndft];
    [s,f,t] = spectrogram(temp_data,g,L,Ndft,SP.Fs);
    figure, waterplot(s,f,t,text),   

    % for M=400:10:800
    %     for L=M-30:M-10
    %         if (L < M)
    %             g = blackman(M); 
    %             %L = 111;    %111
    %             Ndft = max(256,2^nextpow2(M));
    %             texto = [M L Ndft];
    %             [s,f,t] = spectrogram(temp_data,g,L,Ndft,SP.Fs);
    %             waterplot(s,f,t,texto),            
    %             pause(0.001)
    %             % close gcf,
    %         end
    %     end
    % end

    
    figure,
    pspectrum(temp_data,SP.Fs,"spectrogram", ...
    TimeResolution=M/SP.Fs,OverlapPercent=L/M*100, ...
    Leakage=0.7 ...
    ,FrequencyLimits=[0,5e12])
    title("pspectrum")

    %Compute Segment PSDs and Power Spectra
    [s,f,t,p] = spectrogram(temp_data,g,L,Ndft,SP.Fs);
    [r,~,~,q] = spectrogram(temp_data,g,L,Ndft,SP.Fs,"power");
    figure, waterfall(f,t,abs(s)'.^2)
    set(gca,XScale="log",...
    XDir="reverse",View=[30 50])






    
    function waterplot(s,f,t,text)
    % Waterfall plot of spectrogram
        x =waterfall(f,t,abs(s)'.^2);
        set(gca,XDir="reverse",View=[30 50], XLim=[0 5e12] )
        xlabel(text)
        ylabel("Time")
    end
end
