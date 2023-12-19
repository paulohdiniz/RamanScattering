function get_model(SP)
    SP = SP.get_signal_by_time_from_ifft();
    
    meanRealSignal = mean(real(SP.signalIFFT(1).signal),[2 3]);
    meanRealSignal = meanRealSignal(1:length(SP.data_stitched.t_stitched));
    time=SP.data_stitched.t_stitched.';

    %% FUNCTION FIT
    % expr = 'a*cos(b*x + c)*exp(-x/d)';
    % ft = fittype(expr);
    % %3*10^-14*cos(2.5*10^12x + 0.21)*e^(-x*10^11) , in interval (0, 2.5*10^-11)
    % exclude1 =  [time < 8.6e-12 | time > 2.47e-11];
    % fit_result = fit(time, meanRealSignal, ft, ...
    %     'Exclude', exclude1, ...
    %     'StartPoint', [3e-14, 2.5e12, 0.21, 1e-11],...
    %     'TolFun',  1e-25,...
    %     'TolX',  1e-25,...
    %     'DiffMinChange',  1e-25,...,
    %     'MaxFunEvals',  1e25,...,
    %     'MaxIter',  1e3,...
    %     'Display','iter',...
    %     'Algorithm','Levenberg-Marquardt');
    % parametres_sortie = coeffvalues(fit_result)
    % figure, plot(fit_result, time',  squeeze(meanRealSignal)', exclude1);

    %% FUNCTION lsqcurvefit
    time = time';
    meanRealSignal = meanRealSignal';

    %just la descente
    % meanRealSignal= meanRealSignal(time > 8.6e-12 & time < 2.47e-11);
    % time=time(time > 8.6e-12 & time < 2.47e-11);

    fun = @(x, time)x(1)*cos(x(2)*(time - 8.6139e-12) + x(3)).*exp(- (time - 8.6139e-12) * x(4));
    x0 =  [max(meanRealSignal), (12*10^8/SP.c)*2*pi*SP.c, 0.5, SP.peakWidth(1)];
    options = optimoptions('lsqcurvefit',...,
        'OptimalityTolerance', 1e-25, ...
        'FunctionTolerance', 1e-25,...
        'Algorithm','trust-region-reflective',...
        'Display','iter',...,
        'StepTolerance',1e-25,...
        'MaxFunctionEvaluations', 1e3,...
        'MaxIterations', 1e3); %https://stackoverflow.com/questions/45924581/local-minimum-at-initial-point-when-fitting-gaussian-with-lsqcurvefit
    x = lsqcurvefit(fun, x0, time, meanRealSignal, [],[],options);
    figure, plot(time, meanRealSignal, 'ko', time, fun(x0, time), 'b-'), hold on,
    plot(time, meanRealSignal, 'ko', time, fun(x, time), '--')

    %% FUNCTION nlinfit     (Statistics and Machine Learning Toolbox)

    % time = time';
    % meanRealSignal = meanRealSignal';
    % 
    % %just la descente
    % meanRealSignal= meanRealSignal(time > 8.6e-12 & time < 2.47e-11);
    % time=time(time > 8.6e-12 & time < 2.47e-11);
    % 
    % fun = @(x, time)x(1)*cos(x(2)*time + x(3)).*exp(-time / x(4));
    % %x0 =  [3e-13, 2.5e12, 0.5, 7e-10];
    % options = statset('nlinfit');
    % options.TolX =  1e-25;
    % options.TolFun =  1e-25;
    % options.Display = 'iter';
    % options.MaxIter =  1e3;
    % 
    % x = nlinfit(time, meanRealSignal, fun,[], options);
    % figure, plot(time, meanRealSignal, 'ko', time, fun(x0, time), 'b-'), hold on,
    % plot(time, meanRealSignal, 'ko', time, fun(x, time), '--')

    %% Normalized DATA

    % meanRealSignal = meanRealSignal./max(meanRealSignal);
    
    % interval=1:1006;
    % fun = @(x, time)x(1)*cos(x(2)*(time-x(5)) + x(3)).*exp(-abs(time-x(5)) / x(4));
    % x0 =  [5, 100 , 0.5, 0.6, 0.22e-12]; %https://stackoverflow.com/questions/45924581/local-minimum-at-initial-point-when-fitting-gaussian-with-lsqcurvefit
    % options = optimoptions('lsqcurvefit',...,
    %     'OptimalityTolerance', 1e-25, ...
    %     'FunctionTolerance', 1e-25,...
    %     'Algorithm','trust-region-reflective',...
    %     'Display','iter',...,
    %     'StepTolerance',1e-25,...
    %     'MaxFunctionEvaluations', 1e3,...
    %     'MaxIterations', 1e3); %https://stackoverflow.com/questions/45924581/local-minimum-at-initial-point-when-fitting-gaussian-with-lsqcurvefit
    % x = lsqcurvefit(fun, x0, time(interval), meanRealSignal(interval), [],[],options);
    % figure, plot(time, meanRealSignal, 'ko', time, fun(x, time), 'b-')

end

