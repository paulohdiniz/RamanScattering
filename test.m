function test(SP)

    for i=1:3

        filter_window(i).wn =[(SP.peakAmpli_wn(i) - 1.5*(SP.peakWidth(i))) (SP.peakAmpli_wn(i) + 1.5*(SP.peakWidth(i)))];
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
        NFFT = 2^(nextpow2(size(temp,1))); 

        signalIFFT(i).signal = ifft(temp.*SP.hyperspectralRamanImageComplex, NFFT, 1);
        %figure,plot(squeeze(mean(real(signalIFFT(i).signal),[2 3])))

    end
    
    % signal_total = signalIFFT(1).signal + signalIFFT(2).signal + signalIFFT(3).signal;
    % figure, plot(squeeze(mean(real(signal_total),[2 3])));

    meanRealSignal = mean(real(signalIFFT(1).signal),[2 3]);
    meanRealSignal = meanRealSignal(:);
    time = (0:(size(squeeze(meanRealSignal), 1) - 1)) / SP.Fs;
    time = time(:);

    %% Try with fit
    % expr = 'a*cos(b*x + c)*exp(-x/d)';
    % ft = fittype(expr);
    % %3*10^-14*cos(2.5*10^12x + 0.21)*e^(-x*10^11) , in interval (0, 2.5*10^-11)
    % exclude1 =  [time < 8.6e-12 | time > 2.47e-11];
    % fit_result = fit(time, meanRealSignal, ft, ...
    %     'Exclude', exclude1, ...
    %     'StartPoint', [3e-14, 2.5e12, 0.21, 1e-11],...
    %     'TolFun',  1e-20,...
    %     'TolX',  1e-50,...
    %     'DiffMinChange',  1e-20,...,
    %     'MaxFunEvals',  1e20,...,
    %     'MaxIter',  1e20,...
    %     'Display','iter',...
    %     'Algorithm','Levenberg-Marquardt');
    % parametres_sortie = coeffvalues(fit_result)
    % figure, plot(fit_result, time',  squeeze(meanRealSignal)', exclude1);

    %Try with lsqcurvefit
    time = time';
    meanRealSignal = meanRealSignal';

    %just la descente
    meanRealSignal= meanRealSignal(time > 8.6e-12 & time < 2.47e-11);
    time=time(time > 8.6e-12 & time < 2.47e-11);

    fun = @(x, time)x(1)*cos(x(2)*time + x(3)).*exp(-time / x(4));
    x0 =  [3e-13, 2.5e12, 0.5, 7e-10];
    options = optimoptions('lsqcurvefit',...,
        'OptimalityTolerance', 1e-100, ...
        'FunctionTolerance', 1e-100,...
        'Algorithm','trust-region-reflective',...
        'Display','iter',...,
        'StepTolerance',1e-100,...
        'MaxFunctionEvaluations', 1e3,...
        'MaxIterations', 1e4); %https://stackoverflow.com/questions/45924581/local-minimum-at-initial-point-when-fitting-gaussian-with-lsqcurvefit
    x = lsqcurvefit(fun, x0, time, meanRealSignal, [],[],options);
    figure, plot(time, meanRealSignal, 'ko', time, fun(x0, time), 'b-')







end
