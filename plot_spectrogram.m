function plot_spectrogram(SP)

    %Default Values of Spectrogram
    temp_t=SP.data_stitched.t_stitched;
    temp_data=squeeze(mean(SP.data_stitched.data_R,[1 2]));
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
    
    for M=500:10:1000
    % M = 510; %223
    % L= 400;
        for L=M-30:M-1
            if (L < M)
                [M L 2^nextpow2(M)]
                g = hamming(M); 
                %L = 111;    %111
                Ndft = max(1024,2^nextpow2(M));
                [s,f,t] = spectrogram(temp_data,g,L,Ndft,SP.Fs);
                waterplot(s,f,t),            
                pause(0.1)
                % close gcf,
            end
        end

    end

    %Compare spectrogram Function and STFT Definition
    
    % x1 = permute(SP.data_stitched.data_R,[3 1 2]).*repmat(SP.window2.',[1 50 50]);
    % x = squeeze(x1(:, 50, 50));
    % ts = SP.data_stitched.t_stitched;
    % %Use the spectrogram function to compute the STFT of the signal.
    % M = 11;
    % L = 2;
    % g = bartlett(M);
    % Ndft = 2048;
    % 
    % [s,f,t] = spectrogram(x,g,L,Ndft,SP.Fs);
    % waterplot(s,f,t)

    %STFT Definition

    % [segs,~] = buffer(1:length(x),M,L,"nodelay");
    % X = fft(x(segs).*g,Ndft);
    % tbuf = ts(segs);
    % tint = mean(tbuf(2:end,:));
    % fint = 0:SP.Fs/Ndft:SP.Fs-SP.Fs/Ndft;
    % waterplot(X,fint,tint)

    
    function waterplot(s,f,t)
    % Waterfall plot of spectrogram
        x =waterfall(f,t,abs(s)'.^2);
        set(gca,XDir="reverse",View=[30 50], XLim=[0 0.5e13])
        xlabel("Frequency")
        ylabel("Time")
    end
end

