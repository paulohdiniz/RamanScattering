classdef Sid_Processing

    properties

        % Basic Properties and On/Off buttons
        c=3e8;
        DazzlerTimeConversion= 161E-15/1E-6;
        DelayToStartFourier = 1.485e-12;%1.54 for CBZ-DH 1.339E-12
        t0 = 0;
        ReturnOpt = 0;
        ZeroPadOpt = 0;
        stitch = 1;
        DCopt = 1;
        SHGOpt =0;
        Max_wn = 250;

        % Path parameters
        Data_path= pwd;
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
        data_stitched
        t_interp
        data_interp
        window_interp=2

        %FT
        hyperspectralRamanImageComplex
        wn
        ramanSpectrum
        ramanSpectrumWithMask
        rgbimg
        wns_plot
        pixels_plot

        %window
        window2
        window2_name
        ratio_window
        tukey_window_param
        deadtime

        %image processing
        IP = Image_Processing();

    end

    methods

        function Sid_Processing=choose_folders_load_data(Sid_Processing,xp_number)
            Sid_Processing.xp_number = xp_number;
            cd(Sid_Processing.Data_path) % we go there
            folder_content=dir(Sid_Processing.Data_path); % We look a the folder content
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
              
            names=cellfun(@(x) x(9:145),{folder_content(:).name},'UniformOutput',false);
            % % Supprimer la partie finale des éléments de names
            % for i = 1:numel(names)
            %     match = regexp(names{i}, '_([0-9]+)psdelay', 'match');
            %     if ~isempty(match)
            %         names{i} = strrep(names{i}, match, '');
            %     end
            % end
            %TO DO: todos parametros iguais muda o final, para uma serie de exp, o
            %q signifca???
            % for ii=1:length(names)
            %     bla=strcmp(names{i_0},names{ii});
            %     if bla
            %         indice(pp,ii)=ii;
            %     else
            %         pp=pp+1;
            %         i_0=ii;
            %         indice(pp,ii)=ii;
            %     end
            % end
            %Sid_Processing.folders=folder_content(indice(xp_number,indice(xp_number,:)~=0));
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

            end

            % Make the datacubes
            for tt=1:length(Sid_Processing.data_raw)
                Sid_Processing.data_raw(tt).data=reshape(Sid_Processing.data_raw(tt).data,[Sid_Processing.N_t,Sid_Processing.N_x,Sid_Processing.N_y]);
                if Sid_Processing.DCopt
                    Sid_Processing.data_raw(tt).data_R=Sid_Processing.data_raw(tt).data(end:-1:1,1:Sid_Processing.N_x/2,:);
                    Sid_Processing.data_raw(tt).data_T=Sid_Processing.data_raw(tt).data(end:-1:1,Sid_Processing.N_x/2+1:end,:);
                end
            end
        end

        % MAKE A FUNCTION FOR THE WINDOWING OF THE BEGINING PULSE INTERACTION
        function Sid_Processing=stitch_time_axis_no_T(Sid_Processing)
            Total_delay=1/Sid_Processing.Clock_Freq*Sid_Processing.N_t*Sid_Processing.DazzlerTimeConversion;
            time_axis=0:Total_delay/(Sid_Processing.N_t-1):Total_delay;

            for tt=1:length(Sid_Processing.delays)-1
                indice=find(time_axis<(Sid_Processing.delays(tt+1)-Sid_Processing.delays(tt))*1e-12);
                t_delay_parts(tt,:)=Sid_Processing.t0+Sid_Processing.delays(tt)*1e-12+time_axis(indice);
            end
            t_delay_parts(tt+1,:)=Sid_Processing.t0+Sid_Processing.delays(tt+1)*1e-12+time_axis(indice);
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
                %                 if tt==1
                %                      temp_data4D(tt).data=Sid_Processing.data_processed(tt).data_R(1:end-Intersection(tt),:,:);
                %                      temp_data4D(tt).t=time_axis(1:end-Intersection(tt))+Sid_Processing.delays(tt)*1e-12;
                %                 elseif
                %                 tt==length(Sid_Processing.delays)+++
                %                     temp_data4D(tt).data=Sid_Processing.data_processed(tt).data_R(Intersection(tt-1):end,:,:);
                %                      temp_data4D(tt).t=time_axis(Intersection(tt-1):end)+Sid_Processing.delays(tt)*1e-12;
                %                 else
                %                      temp_data4D(tt).data=Sid_Processing.data_processed(tt).data_R(Intersection(tt-1):end-Intersection(tt),:,:);
                %                      temp_data4D(tt).t=time_axis(Intersection(tt-1):end-Intersection(tt))+Sid_Processing.delays(tt)*1e-12;
                %                 end
            end

            t_stitch=[cell2mat({temp_data4D.t})];
            data_R_stitch=permute(cell2mat(cellfun(@(x) permute(x,[2 1 3]),{temp_data4D.data_R},'UniformOutput',false)),[1 3 2]);
            
            % Supprimer les valeurs répétées de t_stitch
            % [unique_t, unique_indices] = unique(t_stitch);
            % data_R_unique = data_R_stitch(:, :, unique_indices);
            % 
            % % Mettre à jour les variables
            % t_stitch = unique_t;
            % data_R_stitch = data_R_unique;
            
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

        function Sid_Processing=window_overlap(Sid_Processing,tukey_window_param,deadtime)
            % We kill the overlap with a half tukey window of the parameter
            % specified
            signal=sum(sum(Sid_Processing.data_raw(1).data_R,2),3);
            N_window=Sid_Processing.N_t*2;
            answer={};
            deadtime2=deadtime;
            while ~any(cellfun(@(x) (strcmp(x,'yes')),answer))
                window=[zeros(deadtime2,1).' tukeywin(N_window,tukey_window_param).' ].';
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
            % We kill the overlap with a half tukey window of the parameter
            % specified
            N_window=Sid_Processing.N_t*2;
            deadtime2=deadtime;
            window=[zeros(deadtime2,1).' tukeywin(N_window,tukey_window_param).' ].';
            window=window(1:Sid_Processing.N_t);

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
        function Sid_Processing=get_values_from_XML(Sid_Processing)
            dom=xmlread('0000_.XML');
            Find = dom.getElementsByTagName('Val');
            Val_ClockFreqValPix=Find.item(75);
            Sid_Processing.Clock_Freq=str2double(char(getNodeValue(Val_ClockFreqValPix.item(0))));
            Val_ScanRangeValPix=Find.item(103);
            Sid_Processing.ScanRange=str2double(char(getNodeValue(Val_ScanRangeValPix.item(0))));

        end
        %         function Sid_Processing=fodlers(Sid_Processing)
        %
        %         end

        %Fourier transform
        function Sid_Processing=FT(Sid_Processing, time_vec, Data_vec)
            dt = diff(time_vec);
            dt = dt(1:end-1);
            dt = mean(abs(dt));
            Fs = 1/dt;
            NFFT = 2^(nextpow2(size(Data_vec,1))+1);
    
            ReducedMeanSignal = bsxfun(@minus,Data_vec,mean(Data_vec,1));
            Y = dt*(fft(ReducedMeanSignal,NFFT,1));
            F = (((0:1/NFFT:1-1/NFFT)*Fs).');
            Sid_Processing.wn = fftshift(F-Fs/2);
            Sid_Processing.wn = Sid_Processing.wn(1:ceil(length(Sid_Processing.wn)/2))/3e8/100; %remove the negatives and put them in the unit cm^-1
            Sid_Processing.hyperspectralRamanImageComplex = Y;
        end

        function Sid_Processing=pick_fourier_window(Sid_Processing, window)
            Sid_Processing.window2_name = window;
            switch window

                %Blackman window
                case 'blackman'
                    Sid_Processing.window2 = [blackman(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                
                %Tukey (tapered cosine) window
                case 'tukeywin'
                    Sid_Processing.window2 = [zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched))) tukeywin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched)),Sid_Processing.tukey_window_param).' ];

                %Hamming window
                case 'hamming'
                    Sid_Processing.window2 = [hamming(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];

                %Hann (Hanning) window
                case 'hann'
                    Sid_Processing.window2 = [hann(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];

                %Flat top weighted window
                case 'flattopwin'
                    Sid_Processing.window2 = [flattopwin(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];

                case '...'
                    Sid_Processing.window2 = [blackman(round(Sid_Processing.ratio_window*length(Sid_Processing.data_stitched.t_stitched))).' zeros(1,round((1-Sid_Processing.ratio_window)*length(Sid_Processing.data_stitched.t_stitched)))];
                    Sid_Processing.window2=ones(size(Sid_Processing.window2));
                   
            end
            Sid_Processing.window2=Sid_Processing.window2(1:length(Sid_Processing.data_stitched.t_stitched));
        end

        function Sid_Processing=make_raman_spectrum(Sid_Processing)
            
            Sid_Processing.ramanSpectrum = abs(sum(sum(Sid_Processing.hyperspectralRamanImageComplex,2),3));
            Sid_Processing.ramanSpectrum = Sid_Processing.ramanSpectrum(1:ceil(length(Sid_Processing.ramanSpectrum)/2));
            
            % %Here it takes the values of the 3 biggest peaks and saves it in the wns_plot variable to plot later
            % [yPeaks,xPeaks,k,prm] = findpeaks(Sid_Processing.ramanSpectrum,Sid_Processing.wn,'SortStr','descend','MinPeakDistance',10);
            % [~,i] = sort(prm,'descend');
            %  if (size(xPeaks) ~= 0)
            %      Sid_Processing.wns_plot=[xPeaks(1) xPeaks(2) xPeaks(3)];
            %  else
            %      Sid_Processing.wns_plot=[0 0 0];
            %  end
            % 
            % 
            % %Here it takes the indices for the plot
            % for i=1:size(Sid_Processing.wns_plot,2)
            %     Sid_Processing.pixels_plot(i) = find(abs(Sid_Processing.wn-Sid_Processing.wns_plot(i))==min(abs(Sid_Processing.wn-Sid_Processing.wns_plot(i))));
            % end
        end

        function Sid_Processing=calculated_ssim_per_wn(Sid_Processing)

            %Reference image
            temp2=squeeze(mean((Sid_Processing.data_processed(1).data_T),[1]));
            temp2=(temp2-min(temp2(:)))./max(temp2(:));
            Sid_Processing.IP.mat_ref = temp2;

            for j = 1:size(Sid_Processing.hyperspectralRamanImageComplex,1)/2
                temp=squeeze(abs(Sid_Processing.hyperspectralRamanImageComplex(j,:,:)));
                temp=(temp-min(temp(:)))./max(temp(:));
                Sid_Processing.IP.mat_img_wn{j} = temp;
                Sid_Processing.IP.ssim_wn(j) = ssim(temp,temp2);
            end

            %TODO : voir pq parfois il arrive ssim de 1024 ele et wn de
            %512.
            if (length(Sid_Processing.IP.ssim_wn) > length(Sid_Processing.wn))
                Sid_Processing.IP.ssim_wn=Sid_Processing.IP.ssim_wn(1:length(Sid_Processing.wn));
            end
            if (length(Sid_Processing.IP.ssim_wn) < length(Sid_Processing.wn))
                Sid_Processing.wn=Sid_Processing.wn(1:length(Sid_Processing.IP.ssim_wn));
            end

        end

        function Sid_Processing=points_to_plot_by_ssim(Sid_Processing)

            temp=Sid_Processing.IP.ssim_wn;
            
            wn_temp=Sid_Processing.wn(Sid_Processing.wn<250);
            temp=(temp(Sid_Processing.wn<250));
            
            wn_temp2=wn_temp(wn_temp>4);
            temp=(temp(wn_temp>4));


            [yPeaks,xPeaks] = findpeaks(temp,wn_temp2,'SortStr','descend','MinPeakHeight',0.05,'MinPeakDistance',2,'MinPeakProminence',0.02);
            
            if numel(xPeaks) ~= 0
                Sid_Processing.wns_plot = zeros(1, numel(xPeaks));
                Sid_Processing.IP.peaks_ssim = zeros(1, numel(yPeaks));
                for i = 1:numel(xPeaks)
                    Sid_Processing.wns_plot(i) = xPeaks(i);
                    Sid_Processing.IP.peaks_ssim(i) = yPeaks(i);
                end
            else
                Sid_Processing.IP.peaks_ssim_wn = [0 0 0];
            end

            %Here it takes the indices for the plot
            for i=1:size(Sid_Processing.wns_plot,2)
                Sid_Processing.pixels_plot(i) = find(abs(Sid_Processing.wn-Sid_Processing.wns_plot(i))==min(abs(Sid_Processing.wn-Sid_Processing.wns_plot(i))));
            end

            if(size(Sid_Processing.wns_plot) ==2)
                Sid_Processing.pixels_plot(3) = Sid_Processing.pixels_plot(2);
            end
                       
            if(size(Sid_Processing.wns_plot) ==1)
                Sid_Processing.pixels_plot(2) = Sid_Processing.pixels_plot(1);
                Sid_Processing.pixels_plot(3) = Sid_Processing.pixels_plot(1);
            end
        end

        function Sid_Processing=points_to_plot_by_frequency(Sid_Processing)

            
            temp=Sid_Processing.ramanSpectrum;
            temp=temp./max(temp(:));
            
            wn_temp=Sid_Processing.wn(Sid_Processing.wn<250);
            temp=(temp(Sid_Processing.wn<250));
            
            wn_temp2=wn_temp(wn_temp>4);
            temp=(temp(wn_temp>4));
            
            %Here it takes the values of the 3 biggest peaks and saves it in the wns_plot variable to plot later
            [yPeaks,xPeaks,k,prm] = findpeaks(temp,wn_temp2,'SortStr','descend','MinPeakHeight',0.05,'MinPeakDistance',2,'MinPeakProminence',0.02);
            [~,i] = sort(prm,'descend');
             if (size(xPeaks) ~= 0)
                 Sid_Processing.wns_plot=[xPeaks(1) xPeaks(2) xPeaks(3)];
             else
                 Sid_Processing.wns_plot=[0 0 0];
             end

            %Here it takes the indices for the plot
            for i=1:size(Sid_Processing.wns_plot,2)
                Sid_Processing.pixels_plot(i) = find(abs(Sid_Processing.wn-Sid_Processing.wns_plot(i))==min(abs(Sid_Processing.wn-Sid_Processing.wns_plot(i))));
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
            newSP.Clock_Freq=obj.Clock_Freq;
            newSP.delays=obj.delays;
            newSP.folders=obj.folders;
            newSP.data_raw=obj.data_raw;
            newSP.data_processed=obj.data_processed;
            newSP.parameters_raw=obj.parameters_raw;
            newSP.N_x=obj.N_x;
            newSP.N_y=obj.N_y;
            newSP.N_t=obj.N_t;
            newSP.data_stitched=obj.data_stitched;
            newSP.t_interp=obj.t_interp;
            newSP.data_interp=obj.data_interp;
            newSP.hyperspectralRamanImageComplex=obj.hyperspectralRamanImageComplex;
            newSP.wn=obj.wn;
            newSP.ramanSpectrum=obj.ramanSpectrum;
            newSP.ramanSpectrumWithMask=obj.ramanSpectrumWithMask;
            newSP.rgbimg=obj.rgbimg;
            newSP.wns_plot=obj.wns_plot;
            newSP.pixels_plot=obj.pixels_plot;
            newSP.window2=obj.window2;
            newSP.window2_name=obj.window2_name;
            newSP.ratio_window=obj.ratio_window;
            newSP.tukey_window_param=obj.tukey_window_param;
            newSP.deadtime=obj.deadtime;
        end
        
    end


end
