function [t_values_for_stop, t_subsections, y_subsections, t, ampli, freq_of_peaks] = get_ideal_times_by_spectrogram(SP)

    %Default Values of Spectrogram
    temp_t=SP.data_stitched.t_stitched;
    temp_data=squeeze(mean(SP.data_stitched.data_R,[1 2]));
    freqLim = 150; %in cm^-1
    percentageOfArea = 0.99;

    M = round(numel(temp_t)*0.6); %pas optimise encore
    L = M - 5 ;    %pas optimise encore, PS: Changer L change seulement la qtd de points de t, s et f ne chageant pas sa taille
    Ndft = max(512,2^nextpow2(M));
    g = blackman(M); 
    [s,f,t] = spectrogram(temp_data,g,L,Ndft,SP.Fs); 

    dt = diff(t);
    freq_of_peaks = arrayfun(@(x) f(find(f >= x * 100 * SP.c, 1)), SP.wns_plot);
    for j=1:numel(freq_of_peaks)
        s_of_peaks{j} = s(f==freq_of_peaks(j),:);
        ampli{j} = abs(s_of_peaks{j})'.^2;

    end

    t_values_for_stop = arrayfun(@(col) t(find(cumtrapz(t, ampli{col})*dt >= percentageOfArea * trapz(t, ampli{col})*dt, 1)), 1:numel(ampli));
    t_subsections = arrayfun(@(col) t(1:find(t == t_values_for_stop(col))), 1:numel(t_values_for_stop), 'UniformOutput', false);
    y_subsections = arrayfun(@(col) ampli{col}(1:find(t == t_values_for_stop(col))), 1:numel(t_values_for_stop), 'UniformOutput', false);
    
    %Si l'amplitude est plus grand qu'1% donc on prendre jusqu'à la fin du temps
    for i=1:numel(ampli)
        if(max(ampli{i})*0.1 <ampli{i}(find(t==t_values_for_stop(i))) )
            t_values_for_stop(i) = temp_t(end); 
        end
    end

end

