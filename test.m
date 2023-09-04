function test(SP)

    for i=1:3

        distance_peak(i).wn =[(SP.peakAmpli_wn(i) - (SP.peakWidth(i)*1.5/2)) (SP.peakAmpli_wn(i) + (SP.peakWidth(i)*1.5/2))]
        
    end
    
    %Here it takes the indices for the plot
    for i=1:size(distance_peak,2)
        inf = find(abs(SP.wn-distance_peak(i).wn(1))==min(abs(SP.wn-distance_peak(i).wn(1))));
        sup = find(abs(SP.wn-distance_peak(i).wn(2))==min(abs(SP.wn-distance_peak(i).wn(2))));
        
        distance_pixels_peak(i).pixels = [inf sup];
    end
    
        hyperSpectralPositive = SP.hyperspectralRamanImageComplex(1:ceil(size(SP.hyperspectralRamanImageComplex, 1) / 2),:,:);
        
        for i =1:size(distance_pixels_peak,2)
            hyperSpectralPeak(i).value = hyperSpectralPositive(distance_pixels_peak(i).pixels(1):distance_pixels_peak(i).pixels(2),:,:);
            signalIFFT = ifft(hyperSpectralPeak(i).value,  size(hyperSpectralPeak(i).value, 1));
            time_signal(i).value =  abs(squeeze(sum(sum(signalIFFT))));
        end
    
    for i = 1:size(distance_pixels_peak, 2)
        figure, hold on,
        time = (0:(size(time_signal(i).value, 1) - 1)) / SP.Fs;
        plot(time, time_signal(i).value);
    end


    model = @(params, t) params(1) * cos(params(2) * t + params(3)) .* exp(-t / params(4));
    %plot 8*10^-11cos(10^14t +  0.7)e^(-t*10^12)H(t) in interval (0, 8*10^-13)
    parametres_initiaux = [8*10e-11, 10e14, 0.7, 10e12];
    parametres_optimaux = lsqcurvefit(model, parametres_initiaux, time, time_signal(2).value');
    A_optimal = parametres_optimaux(1);
    omega_optimal = parametres_optimaux(2);
    phi_optimal = parametres_optimaux(3);
    tau_optimal = parametres_optimaux(4);
    courbe_ajustee = model(parametres_optimaux, time);  % Calculez la courbe ajustée

    plot(time,  time_signal(2).value', '--', time, courbe_ajustee, '-');
    legend('Données', 'Courbe ajustée');
    xlabel('Temps');
    ylabel('Données');
    title('Ajustement du modèle à des données');







end

