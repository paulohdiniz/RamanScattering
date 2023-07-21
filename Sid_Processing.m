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
        % Path parameters
        Data_path='/Users/paulohd/Desktop/PauloDiniz';
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
        % interpolation
        data_stitched
        t_interp
        data_interp
        window_interp=2
    end


    methods

        function Sid_Processing=choose_folders_load_data(Sid_Processing,xp_number)
            cd(Sid_Processing.Data_path) % we go there
            folder_content=dir(Sid_Processing.Data_path); % We look a the folder content
            pp=1;
            i_0=1;
            names=cellfun(@(x) x(9:145),{folder_content(3:end-2).name},'UniformOutput',false); % we avoid the first two.
            for ii=1:length(names)
                bla=strcmp(names{i_0},names{ii});
                if bla
                    indice(pp,ii)=ii;
                else
                    pp=pp+1;
                    i_0=ii;
                    indice(pp,ii)=ii;
                end
            end
            Sid_Processing.folders=folder_content(2+indice(xp_number,indice(xp_number,:)~=0)); % Because we started at 3 before

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
                Sid_Processing.N_t =Sid_Processing.parameters_raw(tt).parameters.data(11);
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
                %                 elseif tt==length(Sid_Processing.delays)
                %                     temp_data4D(tt).data=Sid_Processing.data_processed(tt).data_R(Intersection(tt-1):end,:,:);
                %                      temp_data4D(tt).t=time_axis(Intersection(tt-1):end)+Sid_Processing.delays(tt)*1e-12;
                %                 else
                %                      temp_data4D(tt).data=Sid_Processing.data_processed(tt).data_R(Intersection(tt-1):end-Intersection(tt),:,:);
                %                      temp_data4D(tt).t=time_axis(Intersection(tt-1):end-Intersection(tt))+Sid_Processing.delays(tt)*1e-12;
                %                 end
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

        function Sid_Processing=stitch_interp(Sid_Processing)% to do
            Total_delay=1/SP.Clock_Freq*SP.N_t*SP.DazzlerTimeConversion;
            time_axis=0:Total_delay/(SP.N_t-1):Total_delay;
            t_ini=repmat(time_axis,1,size(SP.delays,2));
            t_ini=t_ini+repelem(SP.delays,264)*1e-12;


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







    end

end
