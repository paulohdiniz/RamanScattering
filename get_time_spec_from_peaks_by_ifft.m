function signalIFFT = get_time_spec_from_peaks_by_ifft(SP)
    

    if (numel(SP.peakAmpli_wn) ~= 0 )
        for i=1:numel(SP.peakAmpli_wn)
    
            filter_window(i).wn =[(SP.peakAmpli_wn(i) - 1.5*(SP.peakWidth(i))) (SP.peakAmpli_wn(i) + 1.5*(SP.peakWidth(i)))];
            complete_filter(i).window = zeros(1, length(SP.ramanSpectrum));
            
        end

    end
    
    for i=1:size(filter_window,2)
        inf = find(abs(SP.wn-filter_window(i).wn(1))==min(abs(SP.wn-filter_window(i).wn(1))));
        sup = find(abs(SP.wn-filter_window(i).wn(2))==min(abs(SP.wn-filter_window(i).wn(2))));
        
        filter_window(i).pixels = [inf sup];
        filter_window(i).window = blackman(sup - inf+1);

        complete_filter(i).window(inf:sup) = ones(1, sup - inf + 1);
        complete_filter(i).window(inf:sup) = complete_filter(i).window(inf:sup)'.*blackman(sup - inf + 1);
        complete_filter(i).window = [complete_filter(i).window complete_filter(i).window(end:-1:1)];

        temp=ones(length(complete_filter(i).window),50,50);
        temp=bsxfun(@times,temp,complete_filter(i).window.');
        NFFT = 2^(nextpow2(size(temp,1))); 

        signalIFFT(i).signal = ifft(temp.*SP.hyperspectralRamanImageComplex, 2048, 1);
        %figure,plot(squeeze(mean(real(signalIFFT(i).signal),[2 3])))

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

