classdef Raman_Processing

    properties

        %Basic Properties and On/Off buttons
        c=299792458;
        DazzlerTimeConversion= 161E-15/1E-6;
        DelayToStartFourier = 1.485e-12;%1.54 for CBZ-DH 1.339E-12
        t0 = 0;
        ReturnOpt = 0;
        ZeroPadOpt = 0;
        stitch = 1;
        DCopt = 1;
        SHGOpt =0;

        %Parameters to find peaks
        Max_wn_to_search_peaks = 250; % < 250cm^-1
        Min_wn_to_search_peaks = 5;   % > 5cm^-1
        MinPeakHeight = 0.05 ; %value normalized, so 5%
        MinPeakDistance = 10 ; %in cm^-1
        MinPeakProminence = 0.02; %value normalized, so 2%

        % Path parameters
        Data_path= cd;
        xp_number;

        % Data parameters
        ScanRange
        Clock_Freq
        delays
        folders
        data_raw
        data_processed
        parameters_raw

        % TF Parameters
        DelayToEndFourier=15E-12;

        %Image Parameters
        N_x
        N_y
        N_t

        %Title Parameters
        list_of_delays
        function_generator
        lockin_parameters
        objectives
        labview_parameters
        laser_parameters
        dazzler_beta2
        polariser_orientation

        % interpolation
        interp_method
        data_stitched
        t_interp
        data_interp
        window_interp=2

        %FT
        hyperspectralRamanImageComplex
        wn
        ramanSpectrum
        ramanSpectrumWithMask

        %Parameters of FFT and spectrogram
        Fs %sampling frequency

        %Properties of Raman Spectrum
        peakAmpli
        peakAmpli_wn
        peakWidth
        peakProm
        scoreCriteria

        %Plot
        wns_plot
        pixels_plot

        %Window
        window2
        window2_name
        ratio_window
        tukey_window_param

        %Impulsion Pulse
        percent_FWHM;
        dead_points;
        
        %image processing
        IP = Image_Processing();

        %PDF
        name_pdf;

        %Model
        signalIFFT;

    end

    methods

        function Raman_Processing=choose_folders_load_data(Raman_Processing,xp_number)
            Raman_Processing.xp_number = xp_number;
            folder_content=dir(pwd);
            folder_content = folder_content(~startsWith({folder_content.name}, ".")); %excludes all starting with "."
            folder_content = folder_content([folder_content.isdir]); %only folders, not files.
            pp=1;
            i_0=1;

            %Taking parameters of title
            full_name_folders=cellfun(@(x) x(1:numel(x)),{folder_content(:).name},'UniformOutput',false);
            for i=1:5
                psdelay = string(regexp(full_name_folders{(xp_number-1)*5 +i}, '_([0-9]+)psdelay', 'tokens')); % match prendre tout
                Raman_Processing.list_of_delays(i).delay = psdelay;
            end
                Raman_Processing.function_generator = string(regexp(full_name_folders{(xp_number-1)*5 +1}, 'FG_([0-9]+)MHZ', 'match', 'once'));
                Raman_Processing.lockin_parameters = string(regexp(full_name_folders{(xp_number-1)*5 +1}, 'APE([\w]+)ns', 'match', 'once'));
                Raman_Processing.objectives = string(regexp(full_name_folders{(xp_number-1)*5 +1}, 'Nik[^_]*(?=_)', 'match', 'once'));
                Raman_Processing.labview_parameters = string(regexp(full_name_folders{(xp_number-1)*5 +1}, '\d*_acc\w*(?=_\d*mW_TiSa)', 'match', 'once'));
                Raman_Processing.laser_parameters = string(regexp(full_name_folders{(xp_number-1)*5 +1}, '\d*mW_TiSa_\d*mW_OPO\d*(?=_)', 'match', 'once'));
                Raman_Processing.dazzler_beta2 = string(regexp(full_name_folders{(xp_number-1)*5 +1}, 'Dazz\w*minus(.*?)(?=_)', 'match', 'once'));
                Raman_Processing.polariser_orientation = string(regexp(full_name_folders{(xp_number-1)*5 +1}, '_([A-Za-z0-9]*QWP\w*?)(?=_([0-9]+)psdelay)', 'tokens', 'once'));
             
                Raman_Processing.folders= folder_content((xp_number-1)*5 +[1 2 3 4 5]);
     
            %Load Data and Parameters
            for tt=1:length(Raman_Processing.folders)
                cd(Raman_Processing.folders(tt).name)
                Raman_Processing.data_raw(tt).data=importdata('0000_.ASC');
                Raman_Processing.parameters_raw(tt).parameters=importdata('0000_.txt');
                Raman_Processing=get_values_from_XML(Raman_Processing);
                cd ..
            end

            % Get the associated Delays
            for tt=1:length(Raman_Processing.folders)
                A=Raman_Processing.folders(tt).name;
                II=strfind(Raman_Processing.folders(tt).name,'_');
                str=erase(A(II(end)+1:end),'psdelay');

                Raman_Processing.delays(tt)=str2double(str);
                Raman_Processing.N_t = Raman_Processing.parameters_raw(tt).parameters.data(11);
                Raman_Processing.N_x = Raman_Processing.parameters_raw(tt).parameters.data(1)*(1+Raman_Processing.DCopt); % In case there is also the transmission data
                Raman_Processing.N_y = Raman_Processing.parameters_raw(tt).parameters.data(2);

                if (isnan(Raman_Processing.delays(tt)))
                    Raman_Processing.delays = [double(struct2array(Raman_Processing.list_of_delays))];
                end

            end

            % Make the datacubes
            for tt=1:length(Raman_Processing.data_raw)
                %Condition to ensure that the reshape will be done correctly
                if(Raman_Processing.N_t*Raman_Processing.N_x*Raman_Processing.N_y ~= size(Raman_Processing.data_raw(tt).data, 1)*size(Raman_Processing.data_raw(tt).data, 2))
                    Raman_Processing.N_y = 50;
                    Raman_Processing.N_x = 100;
                    Raman_Processing.N_t = size(Raman_Processing.data_raw(tt).data, 1)/100;
                end

                Raman_Processing.data_raw(tt).data=reshape(Raman_Processing.data_raw(tt).data,[Raman_Processing.N_t,Raman_Processing.N_x,Raman_Processing.N_y]);
                if Raman_Processing.DCopt
                    Raman_Processing.data_raw(tt).data_R=Raman_Processing.data_raw(tt).data(end:-1:1,1:Raman_Processing.N_x/2,:);
                    Raman_Processing.data_raw(tt).data_T=Raman_Processing.data_raw(tt).data(end:-1:1,Raman_Processing.N_x/2+1:end,:);
                end
            end
        end

        function Raman_Processing=get_values_from_XML(Raman_Processing)
            dom=xmlread('0000_.XML');
            Find = dom.getElementsByTagName('Val');
            Val_ClockFreqValPix=Find.item(75);
            Raman_Processing.Clock_Freq=str2double(char(getNodeValue(Val_ClockFreqValPix.item(0))));
            Val_ScanRangeValPix=Find.item(103);
            Raman_Processing.ScanRange=str2double(char(getNodeValue(Val_ScanRangeValPix.item(0))));

        end

        function Raman_Processing=window_overlap(Raman_Processing,tukey_window_param, pourc_width)
            Total_delay=1/Raman_Processing.Clock_Freq*Raman_Processing.N_t*Raman_Processing.DazzlerTimeConversion;
            time_axis=0:Total_delay/(Raman_Processing.N_t-1):Total_delay;

            first_signal_raw = sum(sum(Raman_Processing.data_raw(1).data_R,2),3);
            signal=sum(sum(Raman_Processing.data_raw(1).data_R,2),3);
            
            first_signal_raw = abs(cumtrapz(first_signal_raw)).^2;
            
            [~,xPeaks, widths, ~] = findpeaks(first_signal_raw,time_axis,'SortStr','descend',...
                            'MinPeakDistance',3e-12);
            
            maxWidthPeak = widths(1);
            time_peak = xPeaks(1);
            answer={};
            while ~any(cellfun(@(x) (strcmp(x,'yes')),answer))

                time_to_zero = time_peak + maxWidthPeak*pourc_width/100;
                dead_points = round((time_to_zero/max(time_axis))*numel(time_axis));
                window=[zeros(dead_points,1).' tukeywin(Raman_Processing.N_t - dead_points,tukey_window_param).' ].'; %TO DO: voir avec sam
    
                window=window(1:Raman_Processing.N_t);
                figure(101),clf,
                plot(signal./max(signal),'DisplayName','Raw signal norm', 'Color', 'blue')
                hold on,plot((signal.*window)./max(signal.*window),'DisplayName','Signal windowed norm','Color', '#006400')
                hold on, plot(window./max(window),"--",'DisplayName','Tukey window norm', 'Color', 'red')
                hold on, plot(first_signal_raw./max(first_signal_raw),'DisplayName','abs(integral of raw)^2', 'Color', 'magenta')
                legend
                prompt = {'Pourcentage of width:','Is the windowing ok ?'};
                dlgtitle = 'Input';
                dims = [1 35];
                definput = {num2str(pourc_width),'no'};
                options.Resize='on';
                options.WindowStyle='normal';
                options.Interpreter='tex';
                answer = inputdlg(prompt,dlgtitle,dims,definput,options);
                pourc_width=str2double(answer{1});
            end
            close(101)
            Raman_Processing.percent_FWHM = pourc_width;

            % THEN AND STITICHING
            for tt=1:size(Raman_Processing.data_raw,2)
                if tt==1
                    Raman_Processing.data_processed(1).data_R=Raman_Processing.data_raw(1).data_R.*repmat(window,1,size(Raman_Processing.data_raw(1).data_R,2),size(Raman_Processing.data_raw(1).data_R,3));
                    Raman_Processing.data_processed(tt).data_T=Raman_Processing.data_raw(tt).data_T;
                else
                    Raman_Processing.data_processed(tt).data_R=Raman_Processing.data_raw(tt).data_R;
                    Raman_Processing.data_processed(tt).data_T=Raman_Processing.data_raw(tt).data_T;
                end
            end
        end

        function Raman_Processing=window_overlap_to_test(Raman_Processing,tukey_window_param, pourc_width) %just for testing

            Total_delay=1/Raman_Processing.Clock_Freq*Raman_Processing.N_t*Raman_Processing.DazzlerTimeConversion;
            time_axis=0:Total_delay/(Raman_Processing.N_t-1):Total_delay;

            first_signal_raw = sum(sum(Raman_Processing.data_raw(1).data_R,2),3);
            
            first_signal_raw = abs(cumtrapz(first_signal_raw)).^2;
            
            [~,xPeaks, widths, ~] = findpeaks(first_signal_raw,time_axis,'SortStr','descend',...
                            'MinPeakDistance',3e-12);
            
            maxWidthPeak = widths(1);
            time_peak = xPeaks(1);

            time_to_zero = time_peak + maxWidthPeak*pourc_width/100;
            Raman_Processing.dead_points = round((time_to_zero/max(time_axis))*numel(time_axis));
            window=[zeros(Raman_Processing.dead_points,1).' tukeywin(Raman_Processing.N_t - Raman_Processing.dead_points,tukey_window_param).' ].'; 
            window=window(1:Raman_Processing.N_t);

            Raman_Processing.percent_FWHM = pourc_width;
            
            % THEN AND STITICHING
            for tt=1:size(Raman_Processing.data_raw,2)
                if tt==1
                    Raman_Processing.data_processed(1).data_R=Raman_Processing.data_raw(1).data_R.*repmat(window,1,size(Raman_Processing.data_raw(1).data_R,2),size(Raman_Processing.data_raw(1).data_R,3));
                    Raman_Processing.data_processed(tt).data_T=Raman_Processing.data_raw(tt).data_T;
                else
                    Raman_Processing.data_processed(tt).data_R=Raman_Processing.data_raw(tt).data_R;
                    Raman_Processing.data_processed(tt).data_T=Raman_Processing.data_raw(tt).data_T;
                end
            end
        end

        function Raman_Processing=Tnorm_and_center_data(Raman_Processing,norm,center) % TO DO

            if norm
                for tt=1:size(Raman_Processing.data_processed,2)
                    Transmission=Raman_Processing.data_processed(tt).data_T;
                    Raman=Raman_Processing.data_processed(tt).data_R;
                    Raman_Processing.data_processed(tt).data_R=Raman./permute(1+permute(Transmission,[2 3 1])./squeeze(max(Transmission,[],1)),[3 1 2]);
                end
            end

            if center
                % The centering only comes for stitching so we do NOT
                % center the first one, just make the others come to the
                % same average pixel by pixel.
                %                 Mean_1=(mean(Raman_Processing.data_processed(1).data_R,1));
                for tt=1:size(Raman_Processing.data_processed,2)
                    mean_per_pixel=repmat(mean(Raman_Processing.data_processed(tt).data_R,1),264,1,1);%-repmat(Mean_1,264,1,1);
                    %                     overall_mean=(mean(Raman_Processing.data_processed(tt).data_R(:)));
                    Raman_Processing.data_processed(tt).data_R=Raman_Processing.data_processed(tt).data_R-mean_per_pixel;
                end
                Raman_Processing.data_processed(1).data_R=Raman_Processing.data_processed(1).data_R-mean(Raman_Processing.data_processed(1).data_R(1:Raman_Processing.dead_points-10));
            end
        end

        function Raman_Processing=stitch_time_axis_T_with_interp(Raman_Processing,interp_method) % TO DO
            Total_delay=1/Raman_Processing.Clock_Freq*Raman_Processing.N_t*Raman_Processing.DazzlerTimeConversion;
            time_axis=0:Total_delay/(Raman_Processing.N_t-1):Total_delay;
            % We find the intersection between the falling first signal and
            % the rising second window.
            Indice_start(1)=1;
            Indice_end(length(Raman_Processing.delays))=Raman_Processing.N_t;
            for tt=1:length(Raman_Processing.delays)-1
                time_delay=(Raman_Processing.delays(tt+1)-Raman_Processing.delays(tt))*1e-12;
                indice_1=find(time_axis>=min(time_axis+time_delay));
                indice_2=find(time_axis+time_delay<=max(time_axis));
                temp1=squeeze(sum(sum(Raman_Processing.data_raw(tt).data_T(indice_1,:,:),2),3));
                temp2=squeeze(sum(sum(Raman_Processing.data_raw(tt+1).data_T(indice_2,:,:),2),3));
                Intersection(tt)=find(abs(temp2-temp1)==min(abs(temp2-temp1))); % THE ERROR IS HERE this is true on the window of 78 so for the beginning it works
                Indice_end(tt)=indice_1(1)+Intersection(tt)-2;
                Indice_start(tt+1)=Intersection(tt);
            end

            % We store them in in a 4D matrix because we don't want to stitch them now.  
            for tt=1:length(Raman_Processing.delays)

                temp_data4D(tt).data_R=Raman_Processing.data_processed(tt).data_R(Indice_start(tt)+Raman_Processing.window_interp:Indice_end(tt)-Raman_Processing.window_interp,:,:);
                temp_data4D(tt).t=time_axis(Indice_start(tt)+Raman_Processing.window_interp:Indice_end(tt)-Raman_Processing.window_interp)+Raman_Processing.delays(tt)*1e-12;
 
            end

            t_stitch=[cell2mat({temp_data4D.t})];
            data_R_stitch=permute(cell2mat(cellfun(@(x) permute(x,[2 1 3]),{temp_data4D.data_R},'UniformOutput',false)),[1 3 2]);
            
            % Interpolate
            t_ini=repmat(time_axis,1,size(Raman_Processing.delays,2));
            t_ini=t_ini+repelem(Raman_Processing.delays,264)*1e-12;
            t_out=0:Total_delay/(Raman_Processing.N_t-1):max(t_ini);

            switch interp_method
                case 'nearest'
                    % ~0.025s interp1 nearest
                    temp= interp1(t_stitch,permute(data_R_stitch, [3 2 1]),t_out,'nearest','extrap');
                    Raman_Processing.data_stitched.data_R = permute(temp,[3 2 1]);
                    Raman_Processing.data_stitched.t_stitched=t_out;

                case 'linear'
                    % ~0.032s interp1 linear
                    temp= interp1(t_stitch,permute(data_R_stitch, [3 2 1]),t_out,'linear','extrap');
                    Raman_Processing.data_stitched.data_R = permute(temp,[3 2 1]);
                    Raman_Processing.data_stitched.t_stitched=t_out;

                case 'makima'
                    temp= makima(t_stitch,data_R_stitch,t_out);
                    Raman_Processing.data_stitched.data_R = temp;
                    Raman_Processing.data_stitched.t_stitched=t_out;

                case 'pchip'
                    temp= pchip(t_stitch,data_R_stitch,t_out);
                    Raman_Processing.data_stitched.data_R = temp;
                    Raman_Processing.data_stitched.t_stitched=t_out;

                case 'spline'
                    temp= spline(t_stitch,data_R_stitch,t_out);
                    Raman_Processing.data_stitched.data_R = temp;
                    Raman_Processing.data_stitched.t_stitched=t_out;
            end
        end

        function Raman_Processing=pick_fourier_window(Raman_Processing, window2_name)
            Raman_Processing.window2_name = window2_name;

            switch window2_name

                %Modified Bartlett-Hann window
                case 'barthannwin'
                    Raman_Processing.window2 = [barthannwin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Bartlett window
                case 'bartlett'
                    Raman_Processing.window2 = [bartlett(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Blackman window
                case 'blackman'
                    Raman_Processing.window2 = [blackman(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                                    
                %Minimum four-term Blackman-Harris window
                case 'blackmanharris'
                    Raman_Processing.window2 = [blackmanharris(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Bohman window
                case 'bohmanwin'
                    Raman_Processing.window2 = [bohmanwin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Chebyshev window
                case 'chebwin'
                    Raman_Processing.window2 = [chebwin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Flat top weighted window
                case 'flattopwin'
                    Raman_Processing.window2 = [flattopwin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Gaussian window
                case 'gausswin'
                    Raman_Processing.window2 = [gausswin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];

                %Hamming window
                case 'hamming'
                    Raman_Processing.window2 = [hamming(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];

                %Hann (Hanning) window
                case 'hann'
                    Raman_Processing.window2 = [hann(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];

                %Kaiser window, we can change beta 
                case 'kaiser'
                    Raman_Processing.window2 = [kaiser(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Nuttall-defined minimum 4-term Blackman-Harris window
                case 'nuttallwin'
                    Raman_Processing.window2 = [nuttallwin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Parzen (de la Vallée Poussin) window
                case 'parzenwin'
                    Raman_Processing.window2 = [parzenwin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Rectangular window
                case 'rectwin'
                    Raman_Processing.window2 = [rectwin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Taylor window
                case 'taylorwin'
                    Raman_Processing.window2 = [taylorwin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
               
                    %Tukey (tapered cosine) window
                case 'tukeywin'
                    Raman_Processing.window2 = [tukeywin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched)),Raman_Processing.tukey_window_param).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                %Tukey (tapered cosine) window
                case 'tukeywinINV'
                    Raman_Processing.window2 = [zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched))) tukeywin(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched)),Raman_Processing.tukey_window_param).' ];
                
                %Triangular window
                case 'triang'
                    Raman_Processing.window2 = [triang(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                
                case 'ones'
                    Raman_Processing.window2 = [blackman(round(Raman_Processing.ratio_window*length(Raman_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Raman_Processing.ratio_window)*length(Raman_Processing.data_stitched.t_stitched)))];
                    Raman_Processing.window2=ones(size(Raman_Processing.window2));
                   
            end
            Raman_Processing.window2=Raman_Processing.window2(1:length(Raman_Processing.data_stitched.t_stitched));
        end

        %Fourier transform
        function Raman_Processing=FT(Raman_Processing, time_vec, Data_vec)
            Raman_Processing=SamplingFrequency(Raman_Processing,time_vec);

            NFFT = 2^(nextpow2(size(Data_vec,1))+1);
            ReducedMeanSignal = bsxfun(@minus,Data_vec,mean(Data_vec,1));
            Y = (1/Raman_Processing.Fs)*(fft(ReducedMeanSignal,NFFT,1));

            F = (((0:1/NFFT:1-1/NFFT)*Raman_Processing.Fs).');
            Raman_Processing.wn = fftshift(F-Raman_Processing.Fs/2);
            Raman_Processing.wn = Raman_Processing.wn(1:ceil(length(Raman_Processing.wn)/2))/Raman_Processing.c/100; %remove the negatives and put them in the unit cm^-1
            Raman_Processing.hyperspectralRamanImageComplex = Y;

        end

        function Raman_Processing=SamplingFrequency(Raman_Processing, time_vec)
            dt = diff(time_vec);
            dt = dt(1:end-1);
            dt = mean(abs(dt));
            Raman_Processing.Fs = 1/dt;
        end

        function Raman_Processing=make_raman_spectrum(Raman_Processing)
            Raman_Processing.ramanSpectrum = abs(sum(sum(Raman_Processing.hyperspectralRamanImageComplex,2),3));
            Raman_Processing.ramanSpectrum = Raman_Processing.ramanSpectrum(1:ceil(length(Raman_Processing.ramanSpectrum)/2));
        end

        function Raman_Processing=calculated_images_scores_per_wn(Raman_Processing)
            %Reference image
            temp2=squeeze(mean((Raman_Processing.data_processed(1).data_T),[1]));
            temp2=(temp2-min(temp2(:)))./max(temp2(:));
            Raman_Processing.IP.mat_ref = temp2;

            for j = 1:size(Raman_Processing.hyperspectralRamanImageComplex,1)/2
                temp=squeeze(abs(Raman_Processing.hyperspectralRamanImageComplex(j,:,:)));
                temp=(temp-min(temp(:)))./max(temp(:));
                Raman_Processing.IP.mat_img_wn{j} = temp;
                Raman_Processing.IP.ssim_wn(j) = ssim(temp,temp2);
                Raman_Processing.IP.piqe_score_wn(j) = piqe(mat2gray(temp));
                Raman_Processing.IP.niqe_score_wn(j) = niqe(mat2gray(temp));
                Raman_Processing.IP.brisque_score_wn(j) = brisque(mat2gray(temp));
            end

            if (length(Raman_Processing.IP.ssim_wn) > length(Raman_Processing.wn))
                Raman_Processing.IP.ssim_wn=Raman_Processing.IP.ssim_wn(1:length(Raman_Processing.wn));
            end
            if (length(Raman_Processing.IP.ssim_wn) < length(Raman_Processing.wn))
                Raman_Processing.wn=Raman_Processing.wn(1:length(Raman_Processing.IP.ssim_wn));
            end

            temp=Raman_Processing.IP.ssim_wn;
            wn_temp=Raman_Processing.wn(Raman_Processing.wn<Raman_Processing.Max_wn_to_search_peaks);
            temp=(temp(Raman_Processing.wn<Raman_Processing.Max_wn_to_search_peaks));  
            wn_temp2=wn_temp(wn_temp>Raman_Processing.Min_wn_to_search_peaks);
            temp=(temp(wn_temp>Raman_Processing.Min_wn_to_search_peaks));

            [yPeaks,xPeaks] = findpeaks(temp,wn_temp2,'SortStr','descend',...
                'MinPeakHeight',Raman_Processing.MinPeakHeight,...
                'MinPeakDistance',Raman_Processing.MinPeakDistance,...
                'MinPeakProminence',Raman_Processing.MinPeakProminence);
            Raman_Processing.IP.peaks_ssim = zeros(1, max(3, numel(yPeaks)));
            Raman_Processing.IP.peaks_ssim_wn = zeros(1, max(3, numel(xPeaks)));
            if numel(xPeaks) ~= 0
                for i = 1:numel(xPeaks)
                    Raman_Processing.IP.peaks_ssim(i) = yPeaks(i);
                    Raman_Processing.IP.peaks_ssim_wn(i) = xPeaks(i);
                end
            else
                Raman_Processing.IP.peaks_ssim_wn = [0 0 0];
                Raman_Processing.IP.peaks_ssim = [0 0 0];
            end

        end

        function Raman_Processing=points_to_plot_by_frequency(Raman_Processing)

            temp=Raman_Processing.ramanSpectrum;
            
            wn_temp=Raman_Processing.wn(Raman_Processing.wn<Raman_Processing.Max_wn_to_search_peaks);
            temp=(temp(Raman_Processing.wn<Raman_Processing.Max_wn_to_search_peaks));
            
            wn_temp2=wn_temp(wn_temp>Raman_Processing.Min_wn_to_search_peaks);
            temp=(temp(wn_temp>Raman_Processing.Min_wn_to_search_peaks));
            
            %Here it takes the values of the 3 biggest peaks and saves it in the wns_plot variable to plot later
            [yPeaks,xPeaks, widths, proms] = findpeaks(temp,wn_temp2,'SortStr','descend',...
                'MinPeakHeight',Raman_Processing.MinPeakHeight*max(temp(:)),...
                'MinPeakDistance',Raman_Processing.MinPeakDistance,...
                'MinPeakProminence',Raman_Processing.MinPeakProminence*max(temp(:)));

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
                Raman_Processing.peakAmpli = zeros(1, numel(yPeaks));
                Raman_Processing.peakAmpli_wn = zeros(1, numel(xPeaks));
                Raman_Processing.peakWidth = zeros(1, numel(xPeaks));
                Raman_Processing.peakProm = zeros(1, numel(xPeaks));
                for i = 1:numel(yPeaks)
                    Raman_Processing.peakAmpli(i) = yPeaksSorted(i);
                    Raman_Processing.peakAmpli_wn(i) = xPeaksSorted(i); %It is in cm^-1
                    Raman_Processing.peakWidth(i) = widthsSorted(i); %It is in cm^-1
                    Raman_Processing.peakProm(i) = promsSorted(i);
                end
                    Raman_Processing.scoreCriteria = get_score_criteria(Raman_Processing);
               
            else
                Raman_Processing.peakAmpli = [0 0 0];
                Raman_Processing.peakAmpli_wn = [0 0 0];
                Raman_Processing.peakWidth = [0 0 0];
                Raman_Processing.peakProm = [0 0 0];
                Raman_Processing.scoreCriteria = [0 0 0 0 0];
            end
             if (size(xPeaks, 1) == 0)
                 Raman_Processing.wns_plot=[0 0 0];
             end
             if (size(xPeaks, 1) == 1)
                 Raman_Processing.wns_plot=[Raman_Processing.peakAmpli_wn(1) 0 0];
             end
             if (size(xPeaks, 1) == 2)
                 Raman_Processing.wns_plot=[Raman_Processing.peakAmpli_wn(1) Raman_Processing.peakAmpli_wn(2) 0];
             end
             if (size(xPeaks, 1) >= 3)
                 Raman_Processing.wns_plot=[Raman_Processing.peakAmpli_wn(1) Raman_Processing.peakAmpli_wn(2) Raman_Processing.peakAmpli_wn(3)];
             end
                 
             
            %Here it takes the indices for the plot
            for i=1:size(Raman_Processing.wns_plot,2)
                Raman_Processing.pixels_plot(i) = find(abs(Raman_Processing.wn-Raman_Processing.wns_plot(i))==min(abs(Raman_Processing.wn-Raman_Processing.wns_plot(i))));
            end
            
            %Case where we have <3 peaks, so we will plot the same image
            if(size(Raman_Processing.wns_plot) ==2)
                Raman_Processing.pixels_plot(3) = Raman_Processing.pixels_plot(2);
                Raman_Processing.wns_plot(3) = Raman_Processing.wns_plot(2);
            end
            if(size(Raman_Processing.wns_plot) ==1)
                Raman_Processing.pixels_plot(2) = Raman_Processing.pixels_plot(1);
                Raman_Processing.pixels_plot(3) = Raman_Processing.pixels_plot(1);
                Raman_Processing.wns_plot(2) = Raman_Processing.pixels_plot(1);
                Raman_Processing.wns_plot(3) = Raman_Processing.pixels_plot(1);
            end
        end

        function Raman_Processing=get_signal_by_time_from_ifft(Raman_Processing)

            if (numel(Raman_Processing.peakAmpli_wn) ~= 0 )
                for i=1:min(numel(Raman_Processing.peakAmpli_wn),3) 
                    filter_window(i).wn =[(Raman_Processing.peakAmpli_wn(i) - 1.5*(Raman_Processing.peakWidth(i))) (Raman_Processing.peakAmpli_wn(i) + 1.5*(Raman_Processing.peakWidth(i)))];
                    complete_filter(i).window = zeros(1, 2*length(Raman_Processing.ramanSpectrum));
                end
            end
            
            for i=1:size(filter_window,2)
                inf = find(abs(Raman_Processing.wn-filter_window(i).wn(1))==min(abs(Raman_Processing.wn-filter_window(i).wn(1))));
                sup = find(abs(Raman_Processing.wn-filter_window(i).wn(2))==min(abs(Raman_Processing.wn-filter_window(i).wn(2))));
                
                filter_window(i).pixels = [inf sup];
                filter_window(i).window = blackman(sup - inf+1);
        
                complete_filter(i).window(inf:sup) = ones(1, sup - inf + 1);
                complete_filter(i).window(inf:sup) = complete_filter(i).window(inf:sup)'.*blackman(sup - inf + 1);
        
                temp=ones(length(complete_filter(i).window),50,50);
                temp=bsxfun(@times,temp,complete_filter(i).window.');
                NFFT = 2^(nextpow2(size(temp,1))); 
                t=(0:(NFFT-1))*(2/Raman_Processing.Fs);
        
                Raman_Processing.signalIFFT(i).signal = ifft(temp.*Raman_Processing.hyperspectralRamanImageComplex(size(complete_filter(i).window(:),1),:,:), NFFT, 1);
        
            end
            if (numel(Raman_Processing.signalIFFT) == 0)
                for i=1:3
                    Raman_Processing.signalIFFT(i).signal = 0;
                end
            end
            if (numel(Raman_Processing.signalIFFT) == 1)
                for i=2:3
                    Raman_Processing.signalIFFT(i).signal = 0;
                end
            end
            if (numel(Raman_Processing.signalIFFT) == 2)
                Raman_Processing.signalIFFT(3).signal = 0;
            end

        end

        function Raman_Processing=make_raman_spectrum_with_mask(Raman_Processing, mask)
        
            % we import a mask 50x50 com 1s and 0s, with the selection made for the user
            % we use the mask to remove these values on hyperspectral
            mask_expanded = repmat(mask, [1, 1, size(Raman_Processing.hyperspectralRamanImageComplex,1)]);
            mask_expanded = permute(mask_expanded, [3, 1, 2]);
            hyperspectralRamanImageComplexFiltered = Raman_Processing.hyperspectralRamanImageComplex.*mask_expanded;

            Raman_Processing.ramanSpectrumWithMask = abs(sum(sum(hyperspectralRamanImageComplexFiltered,2),3));
            Raman_Processing.ramanSpectrumWithMask = Raman_Processing.ramanSpectrumWithMask(1:ceil(length(hyperspectralRamanImageComplexFiltered)/2));
           
        end

        % Deep copy method for RP class        
        function newRP = copy(obj)           
            newRP = Raman_Processing();
            newRP.ScanRange=obj.ScanRange;
            newRP.xp_number=obj.xp_number;
            newRP.Clock_Freq=obj.Clock_Freq;
            newRP.delays=obj.delays;
            newRP.folders=obj.folders;
            newRP.data_raw=obj.data_raw;
            newRP.data_processed=obj.data_processed;
            newRP.parameters_raw=obj.parameters_raw;
            newRP.N_x=obj.N_x;
            newRP.N_y=obj.N_y;
            newRP.N_t=obj.N_t;
            newRP.list_of_delays=obj.list_of_delays;
            newRP.function_generator=obj.function_generator;
            newRP.lockin_parameters=obj.lockin_parameters;
            newRP.objectives=obj.objectives;
            newRP.labview_parameters=obj.labview_parameters;
            newRP.laser_parameters=obj.laser_parameters;
            newRP.dazzler_beta2=obj.dazzler_beta2;
            newRP.polariser_orientation=obj.polariser_orientation;
            newRP.interp_method = obj.interp_method;
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
            newRP.scoreCriteria=obj.scoreCriteria;
            newRP.wns_plot=obj.wns_plot;
            newRP.pixels_plot=obj.pixels_plot;
            newRP.window2=obj.window2;
            newRP.window2_name=obj.window2_name;
            newRP.ratio_window=obj.ratio_window;
            newRP.tukey_window_param=obj.tukey_window_param;
            newRP.percent_FWHM = obj.percent_FWHM;
            newRP.IP = obj.IP;
            newRP.name_pdf = obj.name_pdf;
            newRP.signalIFFT = obj.signalIFFT;
        end
        
    end

end
