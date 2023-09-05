function test(SP)

    for i=1:3

        filter_window(i).wn =[(SP.peakAmpli_wn(i) - (SP.peakWidth(i))) (SP.peakAmpli_wn(i) + (SP.peakWidth(i)))]
        complete_filter(i).window = zeros(1, length(SP.ramanSpectrum));
        
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
        signalIFFT(i).signal = ifft(temp.*SP.hyperspectralRamanImageComplex, 2048, 1);
        %figure,plot(squeeze(mean(real(signalIFFT(i).signal),[2 3])))

    end

    meanRealSignal = mean(real(signalIFFT(1).signal),[2 3]);
    time = (0:(size(squeeze(meanRealSignal), 1) - 1)) / SP.Fs;
    model = @(params, t) params(1) * cos(params(2) * t + params(3)) .* exp(-t / params(4));
    expr = 'a*cos(b*x + c)*exp(-x/d)';
    time = time(:);
    meanRealSignal = meanRealSignal(:);
    exclude1 =  [time < 8.6e-12];
    fit_result = fit(time, meanRealSignal, expr, 'Exclude', exclude1);
    parametres_initiaux = coeffvalues(fit_result);

    lb = []; % Bornes inférieures des paramètres (le cas échéant)
    ub = []; % Bornes supérieures des paramètres (le cas échéant)
    options = optimoptions('lsqcurvefit', 'Display', 'iter'); % Options de l'optimiseur
    fitted_params = lsqcurvefit(@model, parametres_initiaux, time, meanRealSignal, lb, ub, options);

    figure, plot(fit_result, time',  squeeze(meanRealSignal)', exclude1);


    parametres_optimaux = lsqcurvefit(model, parametres_initiaux, time(1:1400), meanRealSignal(1:1400));
    A_optimal = parametres_optimaux(1);
    omega_optimal = parametres_optimaux(2);
    phi_optimal = parametres_optimaux(3);
    tau_optimal = parametres_optimaux(4);
    courbe_ajustee = model(parametres_optimaux, time);  % Calculez la courbe ajustée

    %plot(time,  squeeze(mean(real(signalIFFT(1).signal),[2 3]))', '--');
    plot(time, courbe_ajustee, '-');
    legend('Données', 'Courbe ajustée');
    xlabel('Temps');
    ylabel('Données');
    title('Ajustement du modèle à des données');







end

