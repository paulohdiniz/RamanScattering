function scoreCriteria = get_score_criteria(SP)

       scoreCriteria = [];
       
       %1: peakAmpli.^2
       scoreCriteria(1) = ~isempty(SP.peakAmpli) * sum(SP.peakAmpli(1:min(3, numel(SP.peakAmpli))).^2) + ...
           isempty(SP.peakAmpli) * 0.0;

       %2: peakProm.^2
       scoreCriteria(2) = ~isempty(SP.peakProm) * sum(SP.peakProm(1:min(3, numel(SP.peakProm))).^2) + ...
           isempty(SP.peakProm) * 0.0;

       %2: peakProm./peakWidth
       scoreCriteria(3) = ~isempty(SP.peakProm) * sum(SP.peakProm(1:min(3, numel(SP.peakProm)))./(SP.peakWidth(1:min(3, numel(SP.peakWidth))))) + ...
           isempty(SP.peakProm) * 0.0;
       
       %4: log(peakAmpli) - log(peakWidth)
       scoreCriteria(4) = sum(...
           log(SP.peakAmpli(1:min(3, numel(SP.peakAmpli)))) - ...
           log(SP.peakWidth(1:min(3, numel(SP.peakWidth)))) ...
           ).* ~isempty(SP.peakAmpli) + isempty(SP.peakAmpli) * 0.0;

       %5: log(peakProm) - log(peakWidth)
       scoreCriteria(5) = sum(...
           log(SP.peakProm(1:min(3, numel(SP.peakProm)))) - ...
           log(SP.peakWidth(1:min(3, numel(SP.peakWidth)))) ...
           ).* ~isempty(SP.peakProm) + isempty(SP.peakProm) * 0.0;


end



