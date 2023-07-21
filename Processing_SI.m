classdef Processing_SI

    properties
        %Separate channel/flip time
        Trig_Channel
        Raman_Channel
        Transmission_Channel
        % Trigger_Finder
        indices_trigger
        trig_threshold
        Nt
        shift
        %Filter_param
        filter_order=6
        filter_cutoff
        img_filtered
        % Data_cube
        data_line
        data_cube
        data_windowed
    end

    methods
        function Processing_SI=separate_channels(Processing_SI,img,if_plot)

            switch size(img,3)

                case 1
                    Processing_SI.Transmission_Channel=double(img(:,end:-1:1));
                    if if_plot
                        figure(17),clf,
                        imagesc(squeeze(img))
                        colorbar
                        colormap jet
                        caxis([-2000 2000])
                    end

                case 2
                    Processing_SI.Raman_Channel=double(img(:,end:-1:1,1));
                    Processing_SI.Trig_Channel=double(img(:,end:-1:1,2));

                    if if_plot
                        figure(17),clf,
                        subplot(1,2,1)
                        imagesc(squeeze(img(:,:,1)))
                        colorbar
                        colormap jet
                        caxis([-2000 2000])
                        subplot(1,2,2)
                        imagesc(squeeze(img(:,:,2)))
                        colorbar
                        colormap jet
                        caxis([-2000 2000])

                    end

                case 3
                    Processing_SI.Raman_Channel=double(img(:,end:-1:1,1));
                    Processing_SI.Trig_Channel=double(img(:,end:-1:1,3));
                    Processing_SI.Transmission_Channel=double(img(:,end:-1:1,2));
                    if if_plot
                        figure(17),clf,
                        subplot(1,3,1)
                        imagesc(squeeze(img(:,:,1)))
                        colorbar
                        colormap jet
                        caxis([-2000 2000])
                        subplot(1,3,2)
                        imagesc(squeeze(img(:,:,2)))
                        colorbar
                        colormap jet
                        caxis([-2000 2000])
                        subplot(1,3,3)
                        imagesc(squeeze(img(:,:,3)))
                        colorbar
                        colormap jet
                        caxis([-2000 2000])

                    end

            end
            Processing_SI.data_line=zeros(size(Processing_SI.Raman_Channel,1).^2,Processing_SI.Nt+1);
            Processing_SI.data_cube=zeros(size(Processing_SI.Raman_Channel,1),size(Processing_SI.Raman_Channel,1),Processing_SI.Nt+1);


        end



        function Processing_SI=make_data_cube(Processing_SI)
            temp=Processing_SI.Raman_Channel.';
            temp2=Processing_SI.Trig_Channel.'/10;
            tic
            %  figure(1),clf,plot(temp)
            %  hold on,plot(temp2)
            % temp3=diff(temp2);
            % hold on,plot(1500*(temp3>trig_threshold),'ro')
            Processing_SI=Processing_SI.find_trigger(temp2(:));
            for jj=1:length(Processing_SI.indices_trigger)-1
                Processing_SI.data_line(jj,:)=temp((Processing_SI.indices_trigger(jj)-Processing_SI.shift):(Processing_SI.indices_trigger(jj)+Processing_SI.Nt-Processing_SI.shift));

            end
            Processing_SI.data_line(jj+1,:)=zeros(1,Processing_SI.Nt+1);
            Processing_SI.data_cube=reshape(Processing_SI.data_line,[64 64 Processing_SI.Nt+1]);
            clear Processing_Si.data_line
            toc


        end
        function Processing_SI=filter(Processing_SI)
            [b,a] = butter(Processing_SI.filter_order,Processing_SI.filter_cutoff);
            Processing_SI.img_filtered = filtfilt(b,a,(Processing_SI.Raman_Channel.')).';

        end


        function Processing_SI=find_trigger(Processing_SI,trigger_line)
            trig_diff=diff(trigger_line);
            Processing_SI.indices_trigger=find(trig_diff>Processing_SI.trig_threshold);
            Processing_SI.indices_trigger(Processing_SI.indices_trigger<(Processing_SI.shift+1))=[];
            Processing_SI.indices_trigger(diff(Processing_SI.indices_trigger)<(Processing_SI.Nt-100))=[];
        end


        function Processing_SI=window_overlap(Processing_SI,tukey_window_param,deadtime)
            % We kill the overlap with a half tukey window of the parameter
            % specified
            signal=squeeze(sum(sum(Processing_SI.data_cube,1),2));
            N_window=Processing_SI.Nt*2;

            answer={};
            deadtime2=deadtime;
            while ~any(cellfun(@(x) (strcmp(x,'yes')),answer))
                window=[zeros(deadtime2,1).' tukeywin(N_window,tukey_window_param).' ].';
                window=window(1:Processing_SI.Nt+1);
                figure(101),clf,
                plot(signal)
                hold on,plot(signal.*window,'ro')
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


            Processing_SI.data_windowed=Processing_SI.data_cube.*permute(repmat(window,1,size(Processing_SI.data_cube,1),size(Processing_SI.data_cube,2)),[2 3 1]);
            



        end


    end


end


