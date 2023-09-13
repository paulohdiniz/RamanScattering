classdef Sid_Processing

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
        deadtime

        %image processing
        IP = Image_Processing();

        %PDF
        name_pdf;

    end

    methods

        function Sid_Processing=choose_folders_load_data(Sid_Processing,xp_number)
            Sid_Processing.xp_number = xp_number;
            folder_content=dir(pwd);
            folder_content = folder_content(~startsWith({folder_content.name}, ".")); %excludes all starting with "."
            folder_content = folder_content([folder_content.isdir]); %only folders, not files.
            pp=1;
            i_0=1;

            %Taking parameters of title
            full_name_folders=cellfun(@(x) x(1:numel(x)),{folder_content(:).name},'UniformOutput',false);
            for i=1:5
                psdelay = string(regexp(full_name_folders{(xp_number-1)*5 +i}, '_([0-9]+)psdelay', 'tokens')); % match prendre tout
                Sid_Processing.list_of_delays(i).delay = psdelay;
            end
                Sid_Processing.function_generator = string(regexp(full_name_folders{(xp_number-1)*5 +1}, 'FG_([0-9]+)MHZ', 'match', 'once'));
                Sid_Processing.lockin_parameters = string(regexp(full_name_folders{(xp_number-1)*5 +1}, 'APE([\w]+)ns', 'match', 'once'));
                Sid_Processing.objectives = string(regexp(full_name_folders{(xp_number-1)*5 +1}, 'Nik[^_]*(?=_)', 'match', 'once'));
                Sid_Processing.labview_parameters = string(regexp(full_name_folders{(xp_number-1)*5 +1}, '\d*_acc\w*(?=_\d*mW_TiSa)', 'match', 'once'));
                Sid_Processing.laser_parameters = string(regexp(full_name_folders{(xp_number-1)*5 +1}, '\d*mW_TiSa_\d*mW_OPO\d*(?=_)', 'match', 'once'));
                Sid_Processing.dazzler_beta2 = string(regexp(full_name_folders{(xp_number-1)*5 +1}, 'Dazz\w*minus(.*?)(?=_)', 'match', 'once'));
                Sid_Processing.polariser_orientation = string(regexp(full_name_folders{(xp_number-1)*5 +1}, '_([A-Za-z0-9]*QWP\w*?)(?=_([0-9]+)psdelay)', 'tokens', 'once'));
             
                Sid_Processing.folders= folder_content((xp_number-1)*5 +[1 2 3 4 5]);
     
            %Load Data and Parameters
            for tt=1:length(Sid_Processing.folders)
                cd(Sid_Processing.folders(tt).name)
                Sid_Processing.data_raw(tt).data=importdata('0000_.ASC');
                Sid_Processing.parameters_raw(tt).parameters=importdata('0000_.txt');
                Sid_Processing=get_values_from_XML(Sid_Processing);
                cd ..
            end

            % Get the associated Delays
            for tt=1:length(Sid_Processing.folders)
                A=Sid_Processing.folders(tt).name;
                II=strfind(Sid_Processing.folders(tt).name,'_');
                str=erase(A(II(end)+1:end),'psdelay');

                Sid_Processing.delays(tt)=str2double(str);
                Sid_Processing.N_t = Sid_Processing.parameters_raw(tt).parameters.data(11);
                Sid_Processing.N_x = Sid_Processing.parameters_raw(tt).parameters.data(1)*(1+Sid_Processing.DCopt); % In case there is also the transmission data
                Sid_Processing.N_y = Sid_Processing.parameters_raw(tt).parameters.data(2);

                if (isnan(Sid_Processing.delays(tt)))
                    Sid_Processing.delays = [double(struct2array(Sid_Processing.list_of_delays))];
                end

            end

            % Make the datacubes
            for tt=1:length(Sid_Processing.data_raw)
                %Condition to ensure that the reshape will be done correctly
                if(Sid_Processing.N_t*Sid_Processing.N_x*Sid_Processing.N_y ~= size(Sid_Processing.data_raw(tt).data, 1)*size(Sid_Processing.data_raw(tt).data, 2))
                    Sid_Processing.N_y = 50;
                    Sid_Processing.N_x = 100;
                    Sid_Processing.N_t = size(Sid_Processing.data_raw(tt).data, 1)/100;
                end

                Sid_Processing.data_raw(tt).data=reshape(Sid_Processing.data_raw(tt).data,[Sid_Processing.N_t,Sid_Processing.N_x,Sid_Processing.N_y]);
                if Sid_Processing.DCopt
                    Sid_Processing.data_raw(tt).data_R=Sid_Processing.data_raw(tt).data(end:-1:1,1:Sid_Processing.N_x/2,:);
                    Sid_Processing.data_raw(tt).data_T=Sid_Processing.data_raw(tt).data(end:-1:1,Sid_Processing.N_x/2+1:end,:);
                end
            end
        end

        function Sid_Processing=get_values_from_XML(Sid_Processing)
            dom=xmlread('0000_.XML');
            Find = dom.getElementsByTagName('Val');
            Val_ClockFreqValPix=Find.item(75);
            Sid_Processing.Clock_Freq=str2double(char(getNodeValue(Val_ClockFreqValPix.item(0))));
            Val_ScanRangeValPix=Find.item(103);
            Sid_Processing.ScanRange=str2double(char(getNodeValue(Val_ScanRangeValPix.item(0))));

        end

        function Sid_Processing=window_overlap(Sid_Processing,tukey_window_param,deadtime)
            % We kill the overlap with a half tukey window of the parameter specified
            signal=sum(sum(Sid_Processing.data_raw(1).data_R,2),3);
            answer={};
            deadtime2=deadtime;
            while ~any(cellfun(@(x) (strcmp(x,'yes')),answer))
                window=[zeros(deadtime2,1).' tukeywin(Sid_Processing.N_t - deadtime2,tukey_window_param).' ].'; %TO DO: voir avec sam
                window=window(1:Sid_Processing.N_t); 
                figure(101),clf,
                plot(signal)
                hold on,plot(signal.*window)
                prompt = {'Deadtime:','Is the windowing ok ?'};
                dlgtitle = 'Input';
                dims = [1 35];
                definput = {num2str(deadtime2),'no'};
                options.Resize='on';
                options.WindowStyle='normal';
                options.Interpreter='tex';
                answer = inputdlg(prompt,dlgtitle,dims,definput,options);
                deadtime2=str2double(answer{1});
            end
            close(101)

            % TO DO :
            % THEN AND STITICHING
            for tt=1:size(Sid_Processing.data_raw,2)
                if tt==1
                    Sid_Processing.data_processed(1).data_R=Sid_Processing.data_raw(1).data_R.*repmat(window,1,size(Sid_Processing.data_raw(1).data_R,2),size(Sid_Processing.data_raw(1).data_R,3));
                    Sid_Processing.data_processed(tt).data_T=Sid_Processing.data_raw(tt).data_T;
                else
                    Sid_Processing.data_processed(tt).data_R=Sid_Processing.data_raw(tt).data_R;
                    Sid_Processing.data_processed(tt).data_T=Sid_Processing.data_raw(tt).data_T;
                end
            end
        end

        function Sid_Processing=window_overlap_to_test(Sid_Processing,tukey_window_param,deadtime) %just for testing
            % We kill the overlap with a half tukey window of the parameter specified
            deadtime2=deadtime;
            window=[zeros(deadtime2,1).' tukeywin(Sid_Processing.N_t - deadtime2,tukey_window_param).' ].'; %TO DO: voir avec sam
            window=window(1:Sid_Processing.N_t);

            % THEN AND STITICHING
            for tt=1:size(Sid_Processing.data_raw,2)
                if tt==1
                    Sid_Processing.data_processed(1).data_R=Sid_Processing.data_raw(1).data_R.*repmat(window,1,size(Sid_Processing.data_raw(1).data_R,2),size(Sid_Processing.data_raw(1).data_R,3));
                    Sid_Processing.data_processed(tt).data_T=Sid_Processing.data_raw(tt).data_T;
                else
                    Sid_Processing.data_processed(tt).data_R=Sid_Processing.data_raw(tt).data_R;
                    Sid_Processing.data_processed(tt).data_T=Sid_Processing.data_raw(tt).data_T;
                end
            end
        end

        function Sid_Processing=Tnorm_and_center_data(Sid_Processing,norm,center,deadtime) % TO DO

            if norm
                for tt=1:size(Sid_Processing.data_processed,2)
                    Transmission=Sid_Processing.data_processed(tt).data_T;
                    Raman=Sid_Processing.data_processed(tt).data_R;
                    Sid_Processing.data_processed(tt).data_R=Raman./permute(1+permute(Transmission,[2 3 1])./squeeze(max(Transmission,[],1)),[3 1 2]);
                end
            end

            if center
                % The centering only comes for stitching so we do NOT
                % center the first one, just make the others come to the
                % same average pixel by pixel.
                %                 Mean_1=(mean(Sid_Processing.data_processed(1).data_R,1));
                for tt=1:size(Sid_Processing.data_processed,2)
                    mean_per_pixel=repmat(mean(Sid_Processing.data_processed(tt).data_R,1),264,1,1);%-repmat(Mean_1,264,1,1);
                    %                     overall_mean=(mean(Sid_Processing.data_processed(tt).data_R(:)));
                    Sid_Processing.data_processed(tt).data_R=Sid_Processing.data_processed(tt).data_R-mean_per_pixel;
                end
                Sid_Processing.data_processed(1).data_R=Sid_Processing.data_processed(1).data_R-mean(Sid_Processing.data_processed(1).data_R(1:deadtime-10));
            end
        end

        function Sid_Processing=stitch_time_axis_T_with_interp(Sid_Processing,interp_method) % TO DO
            Total_delay=1/Sid_Processing.Clock_Freq*Sid_Processing.N_t*Sid_Processing.DazzlerTimeConversion;
            time_axis=0:Total_delay/(Sid_Processing.N_t-1):Total_delay;
            % We find the intersection between the falling first signal and
            % the rising second window.
            Indice_start(1)=1;
            Indice_end(length(Sid_Processing.delays))=Sid_Processing.N_t;
            for tt=1:length(Sid_Processing.delays)-1
                time_delay=(Sid_Processing.delays(tt+1)-Sid_Processing.delays(tt))*1e-12;
                indice_1=find(time_axis>=min(time_axis+time_delay));
                indice_2=find(time_axis+time_delay<=max(time_axis));
                temp1=squeeze(sum(sum(Sid_Processing.data_raw(tt).data_T(indice_1,:,:),2),3));
                temp2=squeeze(sum(sum(Sid_Processing.data_raw(tt+1).data_T(indice_2,:,:),2),3));
                Intersection(tt)=find(abs(temp2-temp1)==min(abs(temp2-temp1))); % THE ERROR IS HERE this is true on the window of 78 so for the beginning it works
                Indice_end(tt)=indice_1(1)+Intersection(tt)-2;
                Indice_start(tt+1)=Intersection(tt);
            end

            % We store them in in a 4D matrix because we don't want to stitch them now.  
            for tt=1:length(Sid_Processing.delays)

                temp_data4D(tt).data_R=Sid_Processing.data_processed(tt).data_R(Indice_start(tt)+Sid_Processing.window_interp:Indice_end(tt)-Sid_Processing.window_interp,:,:);
                temp_data4D(tt).t=time_axis(Indice_start(tt)+Sid_Processing.window_interp:Indice_end(tt)-Sid_Processing.window_interp)+Sid_Processing.delays(tt)*1e-12;
 
            end

            t_stitch=[cell2mat({temp_data4D.t})];
            data_R_stitch=permute(cell2mat(cellfun(@(x) permute(x,[2 1 3]),{temp_data4D.data_R},'UniformOutput',false)),[1 3 2]);
            
            % Interpolate
            t_ini=repmat(time_axis,1,size(Sid_Processing.delays,2));
            t_ini=t_ini+repelem(Sid_Processing.delays,264)*1e-12;
            t_out=0:Total_delay/(Sid_Processing.N_t-1):max(t_ini);

            switch interp_method
                case 'makima'
                    temp= makima(t_stitch,data_R_stitch,t_out);
                    Sid_Processing.data_stitched.data_R = temp;
                    Sid_Processing.data_stitched.t_stitched=t_out;

                case 'pchirp'
                    temp= pchip(t_stitch,data_R_stitch,t_out);
                    Sid_Processing.data_stitched.data_R = temp;
                    Sid_Processing.data_stitched.t_stitched=t_out;

                case 'linear'
                    temp= pchip(t_stitch,data_R_stitch,t_out);
                    Sid_Processing.data_stitched.data_R = temp;
                    Sid_Processing.data_stitched.t_stitched=t_out;
            end
        end

        function Sid_Processing=pick_fourier_window(Sid_Processing, window2_name)
            Sid_Processing.window2_name = window2_name;

            switch window2_name

                %Modified Bartlett-Hann window
                case 'barthannwin'
                    Sid_Processing.window2 = [barthannwin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Bartlett window
                case 'bartlett'
                    Sid_Processing.window2 = [bartlett(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Blackman window
                case 'blackman'
                    Sid_Processing.window2 = [blackman(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                                    
                %Minimum four-term Blackman-Harris window
                case 'blackmanharris'
                    Sid_Processing.window2 = [blackmanharris(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Bohman window
                case 'bohmanwin'
                    Sid_Processing.window2 = [bohmanwin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Chebyshev window
                case 'chebwin'
                    Sid_Processing.window2 = [chebwin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Flat top weighted window
                case 'flattopwin'
                    Sid_Processing.window2 = [flattopwin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Gaussian window
                case 'gausswin'
                    Sid_Processing.window2 = [gausswin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];

                %Hamming window
                case 'hamming'
                    Sid_Processing.window2 = [hamming(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];

                %Hann (Hanning) window
                case 'hann'
                    Sid_Processing.window2 = [hann(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];

                %Kaiser window, we can change beta 
                case 'kaiser'
                    Sid_Processing.window2 = [kaiser(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Nuttall-defined minimum 4-term Blackman-Harris window
                case 'nuttallwin'
                    Sid_Processing.window2 = [nuttallwin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Parzen (de la Vallée Poussin) window
                case 'parzenwin'
                    Sid_Processing.window2 = [parzenwin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Rectangular window
                case 'rectwin'
                    Sid_Processing.window2 = [rectwin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Taylor window
                case 'taylorwin'
                    Sid_Processing.window2 = [taylorwin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
               
                    %Tukey (tapered cosine) window
                case 'tukeywin'
                    Sid_Processing.window2 = [tukeywin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched)),Sid_Processing.tukey_window_param).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Tukey (tapered cosine) window
                case 'tukeywinINV'
                    Sid_Processing.window2 = [zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched))) tukeywin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched)),Sid_Processing.tukey_window_param).' ];
                
                %Triangular window
                case 'triang'
                    Sid_Processing.window2 = [triang(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                case 'ones'
                    Sid_Processing.window2 = [blackman(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                    Sid_Processing.window2=ones(size(Sid_Processing.window2));
                   
            end
            Sid_Processing.window2=Sid_Processing.window2(1:length(Sid_Processing.data_stitched.t_stitched));
        end

        %Fourier transform
        function Sid_Processing=FT(Sid_Processing, time_vec, Data_vec)
            Sid_Processing=SamplingFrequency(Sid_Processing,time_vec);

            NFFT = 2^(nextpow2(size(Data_vec,1))+1);
    
            ReducedMeanSignal = bsxfun(@minus,Data_vec,mean(Data_vec,1));
            Y = (1/Sid_Processing.Fs)*(fft(ReducedMeanSignal,NFFT,1));
            varTEMP = fft(ReducedMeanSignal,NFFT,1); %TODO:TEST
            F = (((0:1/NFFT:1-1/NFFT)*Sid_Processing.Fs).');
            Sid_Processing.wn = fftshift(F-Sid_Processing.Fs/2);
            Sid_Processing.wn = Sid_Processing.wn(1:ceil(length(Sid_Processing.wn)/2))/3e8/100; %remove the negatives and put them in the unit cm^-1
            Sid_Processing.hyperspectralRamanImageComplex = Y;
        end

        function Sid_Processing=SamplingFrequency(Sid_Processing, time_vec)
            dt = diff(time_vec);
            dt = dt(1:end-1);
            dt = mean(abs(dt));
            Sid_Processing.Fs = 1/dt;
        end

        function Sid_Processing=make_raman_spectrum(Sid_Processing)
            Sid_Processing.ramanSpectrum = abs(sum(sum(Sid_Processing.hyperspectralRamanImageComplex,2),3));
            Sid_Processing.ramanSpectrum = Sid_Processing.ramanSpectrum(1:ceil(length(Sid_Processing.ramanSpectrum)/2));
        end

        function Sid_Processing=calculated_images_scores_per_wn(Sid_Processing)
            %Reference image
            temp2=squeeze(mean((Sid_Processing.data_processed(1).data_T),[1]));
            temp2=(temp2-min(temp2(:)))./max(temp2(:));
            Sid_Processing.IP.mat_ref = temp2;

            for j = 1:size(Sid_Processing.hyperspectralRamanImageComplex,1)/2
                temp=squeeze(abs(Sid_Processing.hyperspectralRamanImageComplex(j,:,:)));
                temp=(temp-min(temp(:)))./max(temp(:));
                Sid_Processing.IP.mat_img_wn{j} = temp;
                Sid_Processing.IP.ssim_wn(j) = ssim(temp,temp2);
                Sid_Processing.IP.piqe_score_wn(j) = piqe(mat2gray(temp));
                Sid_Processing.IP.niqe_score_wn(j) = niqe(mat2gray(temp));
                Sid_Processing.IP.brisque_score_wn(j) = brisque(mat2gray(temp));
            end

            if (length(Sid_Processing.IP.ssim_wn) > length(Sid_Processing.wn))
                Sid_Processing.IP.ssim_wn=Sid_Processing.IP.ssim_wn(1:length(Sid_Processing.wn));
            end
            if (length(Sid_Processing.IP.ssim_wn) < length(Sid_Processing.wn))
                Sid_Processing.wn=Sid_Processing.wn(1:length(Sid_Processing.IP.ssim_wn));
            end

            temp=Sid_Processing.IP.ssim_wn;
            wn_temp=Sid_Processing.wn(Sid_Processing.wn<Sid_Processing.Max_wn_to_search_peaks);
            temp=(temp(Sid_Processing.wn<Sid_Processing.Max_wn_to_search_peaks));  
            wn_temp2=wn_temp(wn_temp>Sid_Processing.Min_wn_to_search_peaks);
            temp=(temp(wn_temp>Sid_Processing.Min_wn_to_search_peaks));

            [yPeaks,xPeaks] = findpeaks(temp,wn_temp2,'SortStr','descend',...
                'MinPeakHeight',Sid_Processing.MinPeakHeight,...
                'MinPeakDistance',Sid_Processing.MinPeakDistance,...
                'MinPeakProminence',Sid_Processing.MinPeakProminence);
            Sid_Processing.IP.peaks_ssim = zeros(1, max(3, numel(yPeaks)));
            Sid_Processing.IP.peaks_ssim_wn = zeros(1, max(3, numel(xPeaks)));
            if numel(xPeaks) ~= 0
                for i = 1:numel(xPeaks)
                    Sid_Processing.IP.peaks_ssim(i) = yPeaks(i);
                    Sid_Processing.IP.peaks_ssim_wn(i) = xPeaks(i);
                end
            else
                Sid_Processing.IP.peaks_ssim_wn = [0 0 0];
                Sid_Processing.IP.peaks_ssim = [0 0 0];
            end

        end

        function Sid_Processing=points_to_plot_by_ssim(Sid_Processing)

            temp=Sid_Processing.IP.ssim_wn;
            
            wn_temp=Sid_Processing.wn(Sid_Processing.wn<Sid_Processing.Max_wn_to_search_peaks);
            temp=(temp(Sid_Processing.wn<Sid_Processing.Max_wn_to_search_peaks));
            
            wn_temp2=wn_temp(wn_temp>Sid_Processing.Min_wn_to_search_peaks);
            temp=(temp(wn_temp>Sid_Processing.Min_wn_to_search_peaks));


            [yPeaks,xPeaks] = findpeaks(temp,wn_temp2,'SortStr','descend',...
                'MinPeakHeight',Sid_Processing.MinPeakHeight,...
                'MinPeakDistance',Sid_Processing.MinPeakDistance,...
                'MinPeakProminence',Sid_Processing.MinPeakProminence);
            
            if numel(xPeaks) ~= 0
                Sid_Processing.wns_plot = zeros(1, numel(xPeaks));
                for i = 1:numel(xPeaks)
                    Sid_Processing.wns_plot(i) = xPeaks(i);
                end
            else
                Sid_Processing.IP.peaks_ssim_wn = [0 0 0];
            end

            %Here it takes the indices for the plot
            for i=1:size(Sid_Processing.wns_plot,2)
                Sid_Processing.pixels_plot(i) = find(abs(Sid_Processing.wn-Sid_Processing.wns_plot(i))==min(abs(Sid_Processing.wn-Sid_Processing.wns_plot(i))));
            end

            %Case where we have < 3 peaks, so we will plot the same image
            if(size(Sid_Processing.wns_plot) ==2)
                Sid_Processing.pixels_plot(3) = Sid_Processing.pixels_plot(2);
                Sid_Processing.wns_plot(3) = Sid_Processing.wns_plot(2);
            end
            if(size(Sid_Processing.wns_plot) ==1)
                Sid_Processing.pixels_plot(2) = Sid_Processing.pixels_plot(1);
                Sid_Processing.pixels_plot(3) = Sid_Processing.pixels_plot(1);
                Sid_Processing.wns_plot(2) = Sid_Processing.pixels_plot(1);
                Sid_Processing.wns_plot(3) = Sid_Processing.pixels_plot(1);
            end
        end

        function Sid_Processing=points_to_plot_by_frequency(Sid_Processing)

            temp=Sid_Processing.ramanSpectrum;
            temp=temp./max(temp(:));
            
            wn_temp=Sid_Processing.wn(Sid_Processing.wn<Sid_Processing.Max_wn_to_search_peaks);
            temp=(temp(Sid_Processing.wn<Sid_Processing.Max_wn_to_search_peaks));
            
            wn_temp2=wn_temp(wn_temp>Sid_Processing.Min_wn_to_search_peaks);
            temp=(temp(wn_temp>Sid_Processing.Min_wn_to_search_peaks));
            
            %Here it takes the values of the 3 biggest peaks and saves it in the wns_plot variable to plot later
            [yPeaks,xPeaks, widths, proms] = findpeaks(temp,wn_temp2,'SortStr','descend',...
                'MinPeakHeight',Sid_Processing.MinPeakHeight,...
                'MinPeakDistance',Sid_Processing.MinPeakDistance,...
                'MinPeakProminence',Sid_Processing.MinPeakProminence);

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
                Sid_Processing.peakAmpli = zeros(1, numel(yPeaks));
                Sid_Processing.peakAmpli_wn = zeros(1, numel(xPeaks));
                Sid_Processing.peakWidth = zeros(1, numel(xPeaks));
                Sid_Processing.peakProm = zeros(1, numel(xPeaks));
                for i = 1:numel(yPeaks)
                    Sid_Processing.peakAmpli(i) = yPeaksSorted(i);
                    Sid_Processing.peakAmpli_wn(i) = xPeaksSorted(i); %It is in cm^-1
                    Sid_Processing.peakWidth(i) = widthsSorted(i); %It is in cm^-1
                    Sid_Processing.peakProm(i) = promsSorted(i);
                end
                    Sid_Processing.scoreCriteria = get_score_criteria(Sid_Processing);
               
            else
                Sid_Processing.peakAmpli = [0 0 0];
                Sid_Processing.peakAmpli_wn = [0 0 0];
                Sid_Processing.peakWidth = [0 0 0];
                Sid_Processing.peakProm = [0 0 0];
                Sid_Processing.scoreCriteria = [0 0 0 0 0];
            end
             if (size(xPeaks, 1) == 0)
                 Sid_Processing.wns_plot=[0 0 0];
             end
             if (size(xPeaks, 1) == 1)
                 Sid_Processing.wns_plot=[Sid_Processing.peakAmpli_wn(1) 0 0];
             end
             if (size(xPeaks, 1) == 2)
                 Sid_Processing.wns_plot=[Sid_Processing.peakAmpli_wn(1) Sid_Processing.peakAmpli_wn(2) 0];
             end
             if (size(xPeaks, 1) >= 3)
                 Sid_Processing.wns_plot=[Sid_Processing.peakAmpli_wn(1) Sid_Processing.peakAmpli_wn(2) Sid_Processing.peakAmpli_wn(3)];
             end
                 
             
            %Here it takes the indices for the plot
            for i=1:size(Sid_Processing.wns_plot,2)
                Sid_Processing.pixels_plot(i) = find(abs(Sid_Processing.wn-Sid_Processing.wns_plot(i))==min(abs(Sid_Processing.wn-Sid_Processing.wns_plot(i))));
            end
            
            %Case where we have <3 peaks, so we will plot the same image
            if(size(Sid_Processing.wns_plot) ==2)
                Sid_Processing.pixels_plot(3) = Sid_Processing.pixels_plot(2);
                Sid_Processing.wns_plot(3) = Sid_Processing.wns_plot(2);
            end
            if(size(Sid_Processing.wns_plot) ==1)
                Sid_Processing.pixels_plot(2) = Sid_Processing.pixels_plot(1);
                Sid_Processing.pixels_plot(3) = Sid_Processing.pixels_plot(1);
                Sid_Processing.wns_plot(2) = Sid_Processing.pixels_plot(1);
                Sid_Processing.wns_plot(3) = Sid_Processing.pixels_plot(1);
            end
        end

        function Sid_Processing=make_raman_spectrum_with_mask(Sid_Processing, mask)
        
            % we import a mask 50x50 com 1s and 0s, with the selection made for the user
            % we use the mask to remove these values on hyperspectral
            mask_expanded = repmat(mask, [1, 1, size(Sid_Processing.hyperspectralRamanImageComplex,1)]);
            mask_expanded = permute(mask_expanded, [3, 1, 2]);
            hyperspectralRamanImageComplexFiltered = Sid_Processing.hyperspectralRamanImageComplex.*mask_expanded;

            Sid_Processing.ramanSpectrumWithMask = abs(sum(sum(hyperspectralRamanImageComplexFiltered,2),3));
            Sid_Processing.ramanSpectrumWithMask = Sid_Processing.ramanSpectrumWithMask(1:ceil(length(hyperspectralRamanImageComplexFiltered)/2));
           
        end

        % Deep copy method for SP class        
        function newSP = copy(obj)           
            newSP = Sid_Processing();
            newSP.ScanRange=obj.ScanRange;
            newSP.xp_number=obj.xp_number;
            newSP.Clock_Freq=obj.Clock_Freq;
            newSP.delays=obj.delays;
            newSP.folders=obj.folders;
            newSP.data_raw=obj.data_raw;
            newSP.data_processed=obj.data_processed;
            newSP.parameters_raw=obj.parameters_raw;
            newSP.N_x=obj.N_x;
            newSP.N_y=obj.N_y;
            newSP.N_t=obj.N_t;
            newSP.list_of_delays=obj.list_of_delays;
            newSP.function_generator=obj.function_generator;
            newSP.lockin_parameters=obj.lockin_parameters;
            newSP.objectives=obj.objectives;
            newSP.labview_parameters=obj.labview_parameters;
            newSP.laser_parameters=obj.laser_parameters;
            newSP.dazzler_beta2=obj.dazzler_beta2;
            newSP.polariser_orientation=obj.polariser_orientation;
            newSP.interp_method = obj.interp_method;
            newSP.data_stitched=obj.data_stitched;
            newSP.t_interp=obj.t_interp;
            newSP.data_interp=obj.data_interp;
            newSP.hyperspectralRamanImageComplex=obj.hyperspectralRamanImageComplex;
            newSP.wn=obj.wn;
            newSP.ramanSpectrum=obj.ramanSpectrum;
            newSP.ramanSpectrumWithMask=obj.ramanSpectrumWithMask;
            newSP.Fs=obj.Fs;
            newSP.peakAmpli=obj.peakAmpli;
            newSP.peakAmpli_wn=obj.peakAmpli_wn;
            newSP.scoreCriteria=obj.scoreCriteria;
            newSP.wns_plot=obj.wns_plot;
            newSP.pixels_plot=obj.pixels_plot;
            newSP.window2=obj.window2;
            newSP.window2_name=obj.window2_name;
            newSP.ratio_window=obj.ratio_window;
            newSP.tukey_window_param=obj.tukey_window_param;
            newSP.deadtime=obj.deadtime;
            newSP.IP = obj.IP;
            newSP.name_pdf = obj.name_pdf;
        end
        
    end

end
