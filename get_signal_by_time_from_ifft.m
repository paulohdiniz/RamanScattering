function signalIFFT = get_signal_by_time_from_ifft(SP)
    
    t_values_for_stop = get_ideal_times_by_spectrogram(SP);
    if (numel(SP.peakAmpli_wn) ~= 0 )
        for i=1:min(numel(SP.peakAmpli_wn),3) 
            filter_window(i).wn =[(SP.peakAmpli_wn(i) - 1.5*(SP.peakWidth(i))) (SP.peakAmpli_wn(i) + 1.5*(SP.peakWidth(i)))];
            complete_filter(i).window = zeros(1, 2*length(SP.ramanSpectrum));
        end
    end
    
    for i=1:size(filter_window,2)
        inf = find(abs(SP.wn-filter_window(i).wn(1))==min(abs(SP.wn-filter_window(i).wn(1))));
        sup = find(abs(SP.wn-filter_window(i).wn(2))==min(abs(SP.wn-filter_window(i).wn(2))));
        
        filter_window(i).pixels = [inf sup];
        filter_window(i).window = blackman(sup - inf+1);

        complete_filter(i).window(inf:sup) = ones(1, sup - inf + 1);
        complete_filter(i).window(inf:sup) = complete_filter(i).window(inf:sup)'.*blackman(sup - inf + 1);

        temp=ones(length(complete_filter(i).window),50,50);
        temp=bsxfun(@times,temp,complete_filter(i).window.');
        NFFT = 2^(nextpow2(size(temp,1))); 
        t=(0:(NFFT-1))*(2/SP.Fs);

        signalIFFT(i).signal = ifft(temp.*SP.hyperspectralRamanImageComplex(size(complete_filter(i).window(:),1),:,:), NFFT, 1);
        % figure,plot(t(1:1024),real(squeeze(mean(signalIFFT(i).signal(1:1024,:,:),[2 3]))))
        % figure,plot(t,real(squeeze(mean(signalIFFT(i).signal,[2 3]))))


    end
    if (numel(signalIFFT) == 0)
        for i=1:3
            signalIFFT(i).signal = 0;
        end
    end
    if (numel(signalIFFT) == 1)
        for i=2:3
            signalIFFT(i).signal = 0;
        end
    end
    if (numel(signalIFFT) == 2)
        signalIFFT(3).signal = 0;
    end

end

