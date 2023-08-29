function plot_spectrogram(SP)

        Ndft = numel(SP.wn(:));
    
        temp=SP.ramanSpectrum;
        temp=temp./max(temp(:));

        wn_temp=SP.wn(SP.wn<SP.Max_wn_to_search_peaks);
        temp=(temp(SP.wn<SP.Max_wn_to_search_peaks));
        
        wn_temp2=wn_temp(wn_temp>SP.Min_wn_to_search_peaks);
        temp=(temp(wn_temp>SP.Min_wn_to_search_peaks));

        x = chirp(wn_temp2,244,wn_temp2(end),0,"quadratic",0,"convex","complex");
        M = 49;
        L = 11;
        g = bartlett(M);
        [s,f,t] = spectrogram(x,g,L,Ndft,1);
        waterplot(s,f,t)
        function waterplot(s,f,t)
        % Waterfall plot of spectrogram
            waterfall(f,t,abs(s)'.^2)
            set(gca,XDir="reverse",View=[30 50])
            xlabel("Frequency (Hz)")
            ylabel("Time (s)")
        end

        %spectrogram(wn_temp2,'yaxis')


end

