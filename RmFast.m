classdef RmFast %Raman_Fast_Processing

    properties

        % Basic Properties and On/Off buttons
        c=299792458;
        DazzlerTimeConversion= 161E-15/1E-6;
        DelayToStartFourier = 1.485e-12;%1.54 for CBZ-DH 1.339E-12
        DelayToEndFourier=15E-12;

        % Parameters to find peaks
        Max_wn_to_search_peaks = 250; % < 250cm^-1
        Min_wn_to_search_peaks = 5;   % > 5cm^-1
        MinPeakHeight = 0.05 ; %value normalized, so 5%
        MinPeakDistance = 10 ; %in cm^-1
        MinPeakProminence = 0.02; %value normalized, so 2%

        % Data parameters
        Clock_Freq %A RENTRER !
        delays
        data_raw
        data_processed
        time_axis

        % Image Parameters
        N_x
        N_y
        N_t

        % Interpolation
        data_stitched
        t_interp
        data_interp
        window_interp=2

        % FT
        hyperspectralRamanImageComplex
        wn
        ramanSpectrum
        ramanSpectrumWithMask

        % Parameters of FFT and spectrogram
        Fs

        % Properties of Raman Spectrum
        peakAmpli
        peakAmpli_wn
        peakWidth
        peakProm

        % Plot
        wns_plot
        pixels_plot

        % Window
        window2
        window2_name
        ratio_window
        tukey_window_param

        % Impulsion Pulse
        pourc_pulse_width;
        dead_points;
       
    end

    methods

        function RmFast=RmFast(data_raw,N_x,N_y,N_t)
            RmFast.data_raw=data_raw;
            RmFast.N_x=N_x;
            RmFast.N_y=N_y;
            RmFast.N_t=N_t;
            RmFast=RmFast.FT();
        end

        function  RmFast=update_data(RmFast,data_raw,N_x,N_y,N_t)
            RmFast.data_raw=data_raw;
            RmFast.N_x=N_x;
            RmFast.N_y=N_y;
            RmFast.N_t=N_t;
        end

        function RmFast=create_time_axis(RmFast, Clock_Freq)
            RmFast.Clock_Freq = Clock_Freq;
            Total_delay=1/RmFast.Clock_Freq*RmFast.N_t*RmFast.DazzlerTimeConversion;
            RmFast.time_axis=0:Total_delay/(RmFast.N_t-1):Total_delay;
        end

        function RmFast=window_overlap(RmFast,tukey_window_param, pourc_width)

            RmFast.create_time_axis();

            first_signal_raw = sum(sum(RmFast.data_raw(1).data_R,2),3);
            
            first_signal_raw = abs(cumtrapz(first_signal_raw)).^2;
            
            % Fastest method to get the peak width
            [~, maxIndex] = max(first_signal_raw);
            halfMax = 0.5 * first_signal_raw(maxIndex);
            indicesAboveHalfMax = first_signal_raw >= halfMax;
            crossingIndices = find(diff(indicesAboveHalfMax) ~= 0);
            FWHM = RmFast.time_axis(crossingIndices(end)) - RmFast.time_axis(crossingIndices(1));
            time_peak=RmFast.time_axis(maxIndex);
            
            % Old/Slow method to get the peak width
            % [~,xPeaks, widths] = findpeaks(first_signal_raw,RmFast.time_axis,'SortStr','descend',...
            %                 'MinPeakDistance',3e-12);
            % FWHM = widths(1);
            % time_peak = xPeaks(1);

            time_to_zero = time_peak + FWHM*pourc_width/100;
            RmFast.dead_points = round((time_to_zero/max(RmFast.time_axis))*numel(RmFast.time_axis));
            window=[zeros(RmFast.dead_points,1).' tukeywin(RmFast.N_t - RmFast.dead_points,tukey_window_param).' ].'; 
            window=window(1:RmFast.N_t);

            RmFast.pourc_pulse_width = pourc_width;
            
            % THEN AND STITICHING
            for tt=1:size(RmFast.data_raw,2)
                if tt==1
                    RmFast.data_processed(1).data_R=RmFast.data_raw(1).data_R.*repmat(window,1,size(RmFast.data_raw(1).data_R,2),size(RmFast.data_raw(1).data_R,3));
                    RmFast.data_processed(tt).data_T=RmFast.data_raw(tt).data_T;
                else
                    RmFast.data_processed(tt).data_R=RmFast.data_raw(tt).data_R;
                    RmFast.data_processed(tt).data_T=RmFast.data_raw(tt).data_T;
                end
            end
        end

        function RmFast=Tnorm_and_center_data(RmFast)

            for tt=1:size(RmFast.data_processed,2)
                Transmission=RmFast.data_processed(tt).data_T;
                Raman=RmFast.data_processed(tt).data_R;
                RmFast.data_processed(tt).data_R=Raman./permute(1+permute(Transmission,[2 3 1])./squeeze(max(Transmission,[],1)),[3 1 2]);
            end
        end

        function RmFast=stitch_time_axis_T_with_interp(RmFast) % TO DO
            RmFast.create_time_axis();

            % We find the intersection between the falling first signal and
            % the rising second window.
            Indice_start(1)=1;
            Indice_end(length(RmFast.delays))=RmFast.N_t;
            for tt=1:length(RmFast.delays)-1
                time_delay=(RmFast.delays(tt+1)-RmFast.delays(tt))*1e-12;
                indice_1=find(RmFast.time_axis>=min(RmFast.time_axis+time_delay));
                indice_2=find(RmFast.time_axis+time_delay<=max(RmFast.time_axis));
                temp1=squeeze(sum(sum(RmFast.data_raw(tt).data_T(indice_1,:,:),2),3));
                temp2=squeeze(sum(sum(RmFast.data_raw(tt+1).data_T(indice_2,:,:),2),3));
                Intersection(tt)=find(abs(temp2-temp1)==min(abs(temp2-temp1))); % THE ERROR IS HERE this is true on the window of 78 so for the beginning it works
                Indice_end(tt)=indice_1(1)+Intersection(tt)-2;
                Indice_start(tt+1)=Intersection(tt);
            end

            % We store them in in a 4D matrix because we don't want to stitch them now.  
            for tt=1:length(RmFast.delays)

                temp_data4D(tt).data_R=RmFast.data_processed(tt).data_R(Indice_start(tt)+RmFast.window_interp:Indice_end(tt)-RmFast.window_interp,:,:);
                temp_data4D(tt).t=RmFast.time_axis(Indice_start(tt)+RmFast.window_interp:Indice_end(tt)-RmFast.window_interp)+RmFast.delays(tt)*1e-12;
 
            end

            t_stitch=[cell2mat({temp_data4D.t})];
            data_R_stitch=permute(cell2mat(cellfun(@(x) permute(x,[2 1 3]),{temp_data4D.data_R},'UniformOutput',false)),[1 3 2]);
            
            % Interpolate
            t_ini=repmat(RmFast.time_axis,1,size(RmFast.delays,2));
            t_ini=t_ini+repelem(RmFast.delays,264)*1e-12;
            t_out=0:Total_delay/(RmFast.N_t-1):max(t_ini);
            
            % temp= makima(t_stitch,data_R_stitch,t_out);
            % RmFast.data_stitched.data_R = temp;
            % RmFast.data_stitched.t_stitched=t_out;

            %!!!!!
            temp2= interp1(t_stitch,permute(data_R_stitch, [3 1 2]),t_out,'linear','extrap');
            RmFast.data_stitched.data_R = permute(temp2,[3 2 1]);
            RmFast.data_stitched.t_stitched=t_out;

        end

        function RmFast=pick_fourier_window(RmFast, window2_name)
            RmFast.window2_name = window2_name;

            switch window2_name

                %Modified Bartlett-Hann window
                case 'barthannwin'
                    RmFast.window2 = [barthannwin(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                
                %Bartlett window
                case 'bartlett'
                    RmFast.window2 = [bartlett(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                
                %Blackman window
                case 'blackman'
                    RmFast.window2 = [blackman(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                                    
                %Minimum four-term Blackman-Harris window
                case 'blackmanharris'
                    RmFast.window2 = [blackmanharris(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                
                %Bohman window
                case 'bohmanwin'
                    RmFast.window2 = [bohmanwin(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                
                %Chebyshev window
                case 'chebwin'
                    RmFast.window2 = [chebwin(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                
                %Flat top weighted window
                case 'flattopwin'
                    RmFast.window2 = [flattopwin(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                
                %Gaussian window
                case 'gausswin'
                    RmFast.window2 = [gausswin(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];

                %Hamming window
                case 'hamming'
                    RmFast.window2 = [hamming(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];

                %Hann (Hanning) window
                case 'hann'
                    RmFast.window2 = [hann(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];

                %Kaiser window, we can change beta 
                case 'kaiser'
                    RmFast.window2 = [kaiser(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                
                %Nuttall-defined minimum 4-term Blackman-Harris window
                case 'nuttallwin'
                    RmFast.window2 = [nuttallwin(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                
                %Parzen (de la Vallée Poussin) window
                case 'parzenwin'
                    RmFast.window2 = [parzenwin(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                
                %Rectangular window
                case 'rectwin'
                    RmFast.window2 = [rectwin(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                
                %Taylor window
                case 'taylorwin'
                    RmFast.window2 = [taylorwin(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched))).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
               
                    %Tukey (tapered cosine) window
                case 'tukeywin'
                    RmFast.window2 = [tukeywin(round(RmFast.ratio_window*length(RmFast.data_stitched.t_stitched)),RmFast.tukey_window_param).' zeros(1,round((1-RmFast.ratio_window)*length(RmFast.data_stitched.t_stitched)))];
                   
            end
            RmFast.window2=RmFast.window2(1:length(RmFast.data_stitched.t_stitched));
        end

        %Fourier transform
        function RmFast=FT(RmFast, time_vec, Data_vec)
            RmFast=SamplingFrequency(RmFast,time_vec);

            NFFT = 2^(nextpow2(size(Data_vec,1))+1);
            ReducedMeanSignal = bsxfun(@minus,Data_vec,mean(Data_vec,1));
            Y = (1/RmFast.Fs)*(fft(ReducedMeanSignal,NFFT,1));

            F = (((0:1/NFFT:1-1/NFFT)*RmFast.Fs).');
            RmFast.wn = fftshift(F-RmFast.Fs/2);
            RmFast.wn = RmFast.wn(1:ceil(length(RmFast.wn)/2))/RmFast.c/100; %remove the negatives and put them in the unit cm^-1
            RmFast.hyperspectralRamanImageComplex = Y;

        end

        function RmFast=SamplingFrequency(RmFast, time_vec)
            dt = diff(time_vec);
            dt = dt(1:end-1);
            dt = mean(abs(dt));
            RmFast.Fs = 1/dt;
        end

        function RmFast=make_raman_spectrum(RmFast)
            RmFast.ramanSpectrum = abs(sum(sum(RmFast.hyperspectralRamanImageComplex,2),3));
            RmFast.ramanSpectrum = RmFast.ramanSpectrum(1:ceil(length(RmFast.ramanSpectrum)/2));
        end

        function RmFast=points_to_plot_by_frequency(RmFast)

            temp=RmFast.ramanSpectrum;
            
            wn_temp=RmFast.wn(RmFast.wn<RmFast.Max_wn_to_search_peaks);
            temp=(temp(RmFast.wn<RmFast.Max_wn_to_search_peaks));
            
            wn_temp2=wn_temp(wn_temp>RmFast.Min_wn_to_search_peaks);
            temp=(temp(wn_temp>RmFast.Min_wn_to_search_peaks));
            
            %Here it takes the values of the 3 biggest peaks and saves it in the wns_plot variable to plot later
            [yPeaks,xPeaks, widths, proms] = findpeaks(temp,wn_temp2,'SortStr','descend',...
                'MinPeakHeight',RmFast.MinPeakHeight*max(temp(:)),...
                'MinPeakDistance',RmFast.MinPeakDistance,...
                'MinPeakProminence',RmFast.MinPeakProminence*max(temp(:)));

            % Créez une matrice avec les colonnes yPeaks, xPeaks, widths et proms
            peakData = [yPeaks, xPeaks, widths, proms];

            % Triez les données par proms en ordre décroissant
            sortedPeakData = sortrows(peakData, -4); % -4 est l'indice de la colonne des proms

            % Obtenez les valeurs triées dans les variables correspondantes
            yPeaksSorted = sortedPeakData(:, 1);
            xPeaksSorted = sortedPeakData(:, 2);
            widthsSorted = sortedPeakData(:, 3);
            promsSorted = sortedPeakData(:, 4);

            if numel(yPeaks) ~= 0
                RmFast.peakAmpli = zeros(1, numel(yPeaks));
                RmFast.peakAmpli_wn = zeros(1, numel(xPeaks));
                RmFast.peakWidth = zeros(1, numel(xPeaks));
                RmFast.peakProm = zeros(1, numel(xPeaks));
                for i = 1:numel(yPeaks)
                    RmFast.peakAmpli(i) = yPeaksSorted(i);
                    RmFast.peakAmpli_wn(i) = xPeaksSorted(i); %It is in cm^-1
                    RmFast.peakWidth(i) = widthsSorted(i); %It is in cm^-1
                    RmFast.peakProm(i) = promsSorted(i);
                end               
            else
                RmFast.peakAmpli = [0 0 0];
                RmFast.peakAmpli_wn = [0 0 0];
                RmFast.peakWidth = [0 0 0];
                RmFast.peakProm = [0 0 0];
            end
             if (size(xPeaks, 1) == 0)
                 RmFast.wns_plot=[0 0 0];
             end
             if (size(xPeaks, 1) == 1)
                 RmFast.wns_plot=[RmFast.peakAmpli_wn(1) 0 0];
             end
             if (size(xPeaks, 1) == 2)
                 RmFast.wns_plot=[RmFast.peakAmpli_wn(1) RmFast.peakAmpli_wn(2) 0];
             end
             if (size(xPeaks, 1) >= 3)
                 RmFast.wns_plot=[RmFast.peakAmpli_wn(1) RmFast.peakAmpli_wn(2) RmFast.peakAmpli_wn(3)];
             end
                 
            %Here it takes the indices for the plot
            for i=1:size(RmFast.wns_plot,2)
                RmFast.pixels_plot(i) = find(abs(RmFast.wn-RmFast.wns_plot(i))==min(abs(RmFast.wn-RmFast.wns_plot(i))));
            end
            
            %Case where we have <3 peaks, so we will plot the same image
            if(size(RmFast.wns_plot) ==2)
                RmFast.pixels_plot(3) = RmFast.pixels_plot(2);
                RmFast.wns_plot(3) = RmFast.wns_plot(2);
            end
            if(size(RmFast.wns_plot) ==1)
                RmFast.pixels_plot(2) = RmFast.pixels_plot(1);
                RmFast.pixels_plot(3) = RmFast.pixels_plot(1);
                RmFast.wns_plot(2) = RmFast.pixels_plot(1);
                RmFast.wns_plot(3) = RmFast.pixels_plot(1);
            end
        end

        function RmFast=make_raman_spectrum_with_mask(RmFast, mask)
            % we import a mask 50x50 com 1s and 0s, with the selection made for the user
            % we use the mask to remove these values on hyperspectral
            mask_expanded = repmat(mask, [1, 1, size(RmFast.hyperspectralRamanImageComplex,1)]);
            mask_expanded = permute(mask_expanded, [3, 1, 2]);
            hyperspectralRamanImageComplexFiltered = RmFast.hyperspectralRamanImageComplex.*mask_expanded;

            RmFast.ramanSpectrumWithMask = abs(sum(sum(hyperspectralRamanImageComplexFiltered,2),3));
            RmFast.ramanSpectrumWithMask = RmFast.ramanSpectrumWithMask(1:ceil(length(hyperspectralRamanImageComplexFiltered)/2));
           
        end

        % Deep copy method for RP class        
        function newRP = copy(obj)           
            newRP = RmFast();
            newRP.Clock_Freq=obj.Clock_Freq;
            newRP.delays=obj.delays;
            newRP.data_raw=obj.data_raw;
            newRP.data_processed=obj.data_processed;
            newRP.N_x=obj.N_x;
            newRP.N_y=obj.N_y;
            newRP.N_t=obj.N_t;
            newRP.data_stitched=obj.data_stitched;
            newRP.t_interp=obj.t_interp;
            newRP.data_interp=obj.data_interp;
            newRP.hyperspectralRamanImageComplex=obj.hyperspectralRamanImageComplex;
            newRP.wn=obj.wn;
            newRP.ramanSpectrum=obj.ramanSpectrum;
            newRP.ramanSpectrumWithMask=obj.ramanSpectrumWithMask;
            newRP.Fs=obj.Fs;
            newRP.peakAmpli=obj.peakAmpli;
            newRP.peakAmpli_wn=obj.peakAmpli_wn;
            newRP.wns_plot=obj.wns_plot;
            newRP.pixels_plot=obj.pixels_plot;
            newRP.window2=obj.window2;
            newRP.window2_name=obj.window2_name;
            newRP.ratio_window=obj.ratio_window;
            newRP.tukey_window_param=obj.tukey_window_param;
            newRP.pourc_pulse_width = obj.pourc_pulse_width;
        end
        
    end

end
