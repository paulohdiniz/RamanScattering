figure,plot(SP.ramanSpectrum)
peak_pos=SP.pixels_plot(1)
A=blackman(2*peak_pos);
filtre=[A.' zeros(1,length(SP.ramanSpectrum) -length(A))];
hold on,plot(filtre*1.6e-8)

figure,plot(SP.ramanSpectrum.*filtre.')

  % hyperSpectralPositive = SP.hyperspectralRamanImageComplex(1:ceil(size(SP.hyperspectralRamanImageComplex, 1) / 2),:,:);
        
        % for i =1:size(distance_pixels_peak,2)
        filtre2=[filtre filtre(end:-1:1)];
            temp=ones(length(filtre2),50,50);
            temp=bsxfun(@times,temp,filtre2.');
            % hyperSpectralPeak(i).value = hyperSpectralPositive(distance_pixels_peak(i).pixels(1):distance_pixels_peak(i).pixels(2),:,:);
            signalIFFT = ifft(temp.*SP.hyperspectralRamanImageComplex, 2048, 1);
            % time_signal(i).value =  abs(squeeze(sum(sum(signalIFFT))));
        % end


    % signalIFFT = ifft((mean(SP.hyperspectralRamanImageComplex,[2 3])),size(SP.hyperspectralRamanImageComplex,1), 1);
    figure,plot(SP.data_stitched.t_stitched(:)*1E12,signalIFFT(1:1006)*5e13) 
    hold on, plot(SP.data_stitched.t_stitched(:)*1E12,squeeze(mean(SP.data_stitched.data_R,[1 2])), '--b',LineWidth=.5);