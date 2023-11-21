classdef Raman_Fast_Processing

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
        Clock_Freq
        delays
        data_raw
        data_processed

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

        function Raman_Fast_Processing=window_overlap(Raman_Fast_Processing,tukey_window_param, pourc_width)

            Total_delay=1/Raman_Fast_Processing.Clock_Freq*Raman_Fast_Processing.N_t*Raman_Fast_Processing.DazzlerTimeConversion;
            time_axis=0:Total_delay/(Raman_Fast_Processing.N_t-1):Total_delay;

            first_signal_raw = sum(sum(Raman_Fast_Processing.data_raw(1).data_R,2),3);
            
            first_signal_raw = abs(cumtrapz(first_signal_raw)).^2;
            
            [~,xPeaks, widths] = findpeaks(first_signal_raw,time_axis,'SortStr','descend',...
                            'MinPeakDistance',3e-12);
            
            maxWidthPeak = widths(1);
            time_peak = xPeaks(1);

            time_to_zero = time_peak + maxWidthPeak*pourc_width/100;
            Raman_Fast_Processing.dead_points = round((time_to_zero/max(time_axis))*numel(time_axis));
            window=[zeros(Raman_Fast_Processing.dead_points,1).' tukeywin(Raman_Fast_Processing.N_t - Raman_Fast_Processing.dead_points,tukey_window_param).' ].'; 
            window=window(1:Raman_Fast_Processing.N_t);

            Raman_Fast_Processing.pourc_pulse_width = pourc_width;
            
            % THEN AND STITICHING
            for tt=1:size(Raman_Fast_Processing.data_raw,2)
                if tt==1
                    Raman_Fast_Processing.data_processed(1).data_R=Raman_Fast_Processing.data_raw(1).data_R.*repmat(window,1,size(Raman_Fast_Processing.data_raw(1).data_R,2),size(Raman_Fast_Processing.data_raw(1).data_R,3));
                    Raman_Fast_Processing.data_processed(tt).data_T=Raman_Fast_Processing.data_raw(tt).data_T;
                else
                    Raman_Fast_Processing.data_processed(tt).data_R=Raman_Fast_Processing.data_raw(tt).data_R;
                    Raman_Fast_Processing.data_processed(tt).data_T=Raman_Fast_Processing.data_raw(tt).data_T;
                end
            end
        end

        function Raman_Fast_Processing=Tnorm_and_center_data(Raman_Fast_Processing)

            for tt=1:size(Raman_Fast_Processing.data_processed,2)
                Transmission=Raman_Fast_Processing.data_processed(tt).data_T;
                Raman=Raman_Fast_Processing.data_processed(tt).data_R;
                Raman_Fast_Processing.data_processed(tt).data_R=Raman./permute(1+permute(Transmission,[2 3 1])./squeeze(max(Transmission,[],1)),[3 1 2]);
            end
        end

        function Raman_Fast_Processing=stitch_time_axis_T_with_interp(Raman_Fast_Processing) % TO DO
            Total_delay=1/Raman_Fast_Processing.Clock_Freq*Raman_Fast_Processing.N_t*Raman_Fast_Processing.DazzlerTimeConversion;
            time_axis=0:Total_delay/(Raman_Fast_Processing.N_t-1):Total_delay;
            % We find the intersection between the falling first signal and
            % the rising second window.
            Indice_start(1)=1;
            Indice_end(length(Raman_Fast_Processing.delays))=Raman_Fast_Processing.N_t;
            for tt=1:length(Raman_Fast_Processing.delays)-1
                time_delay=(Raman_Fast_Processing.delays(tt+1)-Raman_Fast_Processing.delays(tt))*1e-12;
                indice_1=find(time_axis>=min(time_axis+time_delay));
                indice_2=find(time_axis+time_delay<=max(time_axis));
                temp1=squeeze(sum(sum(Raman_Fast_Processing.data_raw(tt).data_T(indice_1,:,:),2),3));
                temp2=squeeze(sum(sum(Raman_Fast_Processing.data_raw(tt+1).data_T(indice_2,:,:),2),3));
                Intersection(tt)=find(abs(temp2-temp1)==min(abs(temp2-temp1))); % THE ERROR IS HERE this is true on the window of 78 so for the beginning it works
                Indice_end(tt)=indice_1(1)+Intersection(tt)-2;
                Indice_start(tt+1)=Intersection(tt);
            end

            % We store them in in a 4D matrix because we don't want to stitch them now.  
            for tt=1:length(Raman_Fast_Processing.delays)

                temp_data4D(tt).data_R=Raman_Fast_Processing.data_processed(tt).data_R(Indice_start(tt)+Raman_Fast_Processing.window_interp:Indice_end(tt)-Raman_Fast_Processing.window_interp,:,:);
                temp_data4D(tt).t=time_axis(Indice_start(tt)+Raman_Fast_Processing.window_interp:Indice_end(tt)-Raman_Fast_Processing.window_interp)+Raman_Fast_Processing.delays(tt)*1e-12;
 
            end

            t_stitch=[cell2mat({temp_data4D.t})];
            data_R_stitch=permute(cell2mat(cellfun(@(x) permute(x,[2 1 3]),{temp_data4D.data_R},'UniformOutput',false)),[1 3 2]);
            
            % Interpolate
            t_ini=repmat(time_axis,1,size(Raman_Fast_Processing.delays,2));
            t_ini=t_ini+repelem(Raman_Fast_Processing.delays,264)*1e-12;
            t_out=0:Total_delay/(Raman_Fast_Processing.N_t-1):max(t_ini);
            
            % temp= makima(t_stitch,data_R_stitch,t_out);
            % Raman_Fast_Processing.data_stitched.data_R = temp;
            % Raman_Fast_Processing.data_stitched.t_stitched=t_out;

            %TODO: MAKE interpolation LINEAR !!!!!
            temp2= interp1(t_stitch,permute(data_R_stitch, [3 1 2]),t_out,'linear','extrap');
            Raman_Fast_Processing.data_stitched.data_R = permute(temp2,[3 2 1]);
            Raman_Fast_Processing.data_stitched.t_stitched=t_out;

        end

        function Raman_Fast_Processing=pick_fourier_window(Raman_Fast_Processing, window2_name)
            Raman_Fast_Processing.window2_name = window2_name;

            switch window2_name

                %Modified Bartlett-Hann window
                case 'barthannwin'
                    Raman_Fast_Processing.window2 = [barthannwin(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                
                %Bartlett window
                case 'bartlett'
                    Raman_Fast_Processing.window2 = [bartlett(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                
                %Blackman window
                case 'blackman'
                    Raman_Fast_Processing.window2 = [blackman(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                                    
                %Minimum four-term Blackman-Harris window
                case 'blackmanharris'
                    Raman_Fast_Processing.window2 = [blackmanharris(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                
                %Bohman window
                case 'bohmanwin'
                    Raman_Fast_Processing.window2 = [bohmanwin(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                
                %Chebyshev window
                case 'chebwin'
                    Raman_Fast_Processing.window2 = [chebwin(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                
                %Flat top weighted window
                case 'flattopwin'
                    Raman_Fast_Processing.window2 = [flattopwin(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                
                %Gaussian window
                case 'gausswin'
                    Raman_Fast_Processing.window2 = [gausswin(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];

                %Hamming window
                case 'hamming'
                    Raman_Fast_Processing.window2 = [hamming(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];

                %Hann (Hanning) window
                case 'hann'
                    Raman_Fast_Processing.window2 = [hann(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];

                %Kaiser window, we can change beta 
                case 'kaiser'
                    Raman_Fast_Processing.window2 = [kaiser(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                
                %Nuttall-defined minimum 4-term Blackman-Harris window
                case 'nuttallwin'
                    Raman_Fast_Processing.window2 = [nuttallwin(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                
                %Parzen (de la Vallée Poussin) window
                case 'parzenwin'
                    Raman_Fast_Processing.window2 = [parzenwin(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                
                %Rectangular window
                case 'rectwin'
                    Raman_Fast_Processing.window2 = [rectwin(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                
                %Taylor window
                case 'taylorwin'
                    Raman_Fast_Processing.window2 = [taylorwin(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
               
                    %Tukey (tapered cosine) window
                case 'tukeywin'
                    Raman_Fast_Processing.window2 = [tukeywin(round(Raman_Fast_Processing.ratio_window*length(Raman_Fast_Processing.data_stitched.t_stitched)),Raman_Fast_Processing.tukey_window_param).' zeros(1,round((1-Raman_Fast_Processing.ratio_window)*length(Raman_Fast_Processing.data_stitched.t_stitched)))];
                   
            end
            Raman_Fast_Processing.window2=Raman_Fast_Processing.window2(1:length(Raman_Fast_Processing.data_stitched.t_stitched));
        end

        %Fourier transform
        function Raman_Fast_Processing=FT(Raman_Fast_Processing, time_vec, Data_vec)
            Raman_Fast_Processing=SamplingFrequency(Raman_Fast_Processing,time_vec);

            NFFT = 2^(nextpow2(size(Data_vec,1))+1);
            ReducedMeanSignal = bsxfun(@minus,Data_vec,mean(Data_vec,1));
            Y = (1/Raman_Fast_Processing.Fs)*(fft(ReducedMeanSignal,NFFT,1));

            F = (((0:1/NFFT:1-1/NFFT)*Raman_Fast_Processing.Fs).');
            Raman_Fast_Processing.wn = fftshift(F-Raman_Fast_Processing.Fs/2);
            Raman_Fast_Processing.wn = Raman_Fast_Processing.wn(1:ceil(length(Raman_Fast_Processing.wn)/2))/Raman_Fast_Processing.c/100; %remove the negatives and put them in the unit cm^-1
            Raman_Fast_Processing.hyperspectralRamanImageComplex = Y;

        end

        function Raman_Fast_Processing=SamplingFrequency(Raman_Fast_Processing, time_vec)
            dt = diff(time_vec);
            dt = dt(1:end-1);
            dt = mean(abs(dt));
            Raman_Fast_Processing.Fs = 1/dt;
        end

        function Raman_Fast_Processing=make_raman_spectrum(Raman_Fast_Processing)
            Raman_Fast_Processing.ramanSpectrum = abs(sum(sum(Raman_Fast_Processing.hyperspectralRamanImageComplex,2),3));
            Raman_Fast_Processing.ramanSpectrum = Raman_Fast_Processing.ramanSpectrum(1:ceil(length(Raman_Fast_Processing.ramanSpectrum)/2));
        end

        function Raman_Fast_Processing=points_to_plot_by_frequency(Raman_Fast_Processing)

            temp=Raman_Fast_Processing.ramanSpectrum;
            
            wn_temp=Raman_Fast_Processing.wn(Raman_Fast_Processing.wn<Raman_Fast_Processing.Max_wn_to_search_peaks);
            temp=(temp(Raman_Fast_Processing.wn<Raman_Fast_Processing.Max_wn_to_search_peaks));
            
            wn_temp2=wn_temp(wn_temp>Raman_Fast_Processing.Min_wn_to_search_peaks);
            temp=(temp(wn_temp>Raman_Fast_Processing.Min_wn_to_search_peaks));
            
            %Here it takes the values of the 3 biggest peaks and saves it in the wns_plot variable to plot later
            [yPeaks,xPeaks, widths, proms] = findpeaks(temp,wn_temp2,'SortStr','descend',...
                'MinPeakHeight',Raman_Fast_Processing.MinPeakHeight*max(temp(:)),...
                'MinPeakDistance',Raman_Fast_Processing.MinPeakDistance,...
                'MinPeakProminence',Raman_Fast_Processing.MinPeakProminence*max(temp(:)));

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
                Raman_Fast_Processing.peakAmpli = zeros(1, numel(yPeaks));
                Raman_Fast_Processing.peakAmpli_wn = zeros(1, numel(xPeaks));
                Raman_Fast_Processing.peakWidth = zeros(1, numel(xPeaks));
                Raman_Fast_Processing.peakProm = zeros(1, numel(xPeaks));
                for i = 1:numel(yPeaks)
                    Raman_Fast_Processing.peakAmpli(i) = yPeaksSorted(i);
                    Raman_Fast_Processing.peakAmpli_wn(i) = xPeaksSorted(i); %It is in cm^-1
                    Raman_Fast_Processing.peakWidth(i) = widthsSorted(i); %It is in cm^-1
                    Raman_Fast_Processing.peakProm(i) = promsSorted(i);
                end               
            else
                Raman_Fast_Processing.peakAmpli = [0 0 0];
                Raman_Fast_Processing.peakAmpli_wn = [0 0 0];
                Raman_Fast_Processing.peakWidth = [0 0 0];
                Raman_Fast_Processing.peakProm = [0 0 0];
            end
             if (size(xPeaks, 1) == 0)
                 Raman_Fast_Processing.wns_plot=[0 0 0];
             end
             if (size(xPeaks, 1) == 1)
                 Raman_Fast_Processing.wns_plot=[Raman_Fast_Processing.peakAmpli_wn(1) 0 0];
             end
             if (size(xPeaks, 1) == 2)
                 Raman_Fast_Processing.wns_plot=[Raman_Fast_Processing.peakAmpli_wn(1) Raman_Fast_Processing.peakAmpli_wn(2) 0];
             end
             if (size(xPeaks, 1) >= 3)
                 Raman_Fast_Processing.wns_plot=[Raman_Fast_Processing.peakAmpli_wn(1) Raman_Fast_Processing.peakAmpli_wn(2) Raman_Fast_Processing.peakAmpli_wn(3)];
             end
                 
            %Here it takes the indices for the plot
            for i=1:size(Raman_Fast_Processing.wns_plot,2)
                Raman_Fast_Processing.pixels_plot(i) = find(abs(Raman_Fast_Processing.wn-Raman_Fast_Processing.wns_plot(i))==min(abs(Raman_Fast_Processing.wn-Raman_Fast_Processing.wns_plot(i))));
            end
            
            %Case where we have <3 peaks, so we will plot the same image
            if(size(Raman_Fast_Processing.wns_plot) ==2)
                Raman_Fast_Processing.pixels_plot(3) = Raman_Fast_Processing.pixels_plot(2);
                Raman_Fast_Processing.wns_plot(3) = Raman_Fast_Processing.wns_plot(2);
            end
            if(size(Raman_Fast_Processing.wns_plot) ==1)
                Raman_Fast_Processing.pixels_plot(2) = Raman_Fast_Processing.pixels_plot(1);
                Raman_Fast_Processing.pixels_plot(3) = Raman_Fast_Processing.pixels_plot(1);
                Raman_Fast_Processing.wns_plot(2) = Raman_Fast_Processing.pixels_plot(1);
                Raman_Fast_Processing.wns_plot(3) = Raman_Fast_Processing.pixels_plot(1);
            end
        end

        function Raman_Fast_Processing=make_raman_spectrum_with_mask(Raman_Fast_Processing, mask)
            % we import a mask 50x50 com 1s and 0s, with the selection made for the user
            % we use the mask to remove these values on hyperspectral
            mask_expanded = repmat(mask, [1, 1, size(Raman_Fast_Processing.hyperspectralRamanImageComplex,1)]);
            mask_expanded = permute(mask_expanded, [3, 1, 2]);
            hyperspectralRamanImageComplexFiltered = Raman_Fast_Processing.hyperspectralRamanImageComplex.*mask_expanded;

            Raman_Fast_Processing.ramanSpectrumWithMask = abs(sum(sum(hyperspectralRamanImageComplexFiltered,2),3));
            Raman_Fast_Processing.ramanSpectrumWithMask = Raman_Fast_Processing.ramanSpectrumWithMask(1:ceil(length(hyperspectralRamanImageComplexFiltered)/2));
           
        end

        % Deep copy method for RP class        
        function newRP = copy(obj)           
            newRP = Raman_Fast_Processing();
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
