close all
clear variables 
cd('/Users/paulohd/Desktop/PauloDiniz')
%addpath(genpath('/Users/paulohd/Desktop/PauloDiniz'))
addpath('/Users/paulohd/Desktop/PauloDiniz')
%% CODE FOR USING THE CLASS THAT MAKES THE SAME THING AS SID'S CODE
% We generate the class
SP=Sid_Processing();
SP.DCopt=1;

% Choose the experiment, and load all the data from the correct folders
SP=SP.choose_folders_load_data(3); %1,2,3


SP.ratio_window = 1;
SP.tukey_window_param = 0.3;
SP.deadtime=79;

SP=SP.window_overlap(SP.tukey_window_param,SP.deadtime);

SP=SP.Tnorm_and_center_data(1,0,SP.deadtime);

SP=SP.stitch_time_axis_T_with_interp('makima'); %makima, pchirp, linear

SP = SP.pick_window('...'); %blackman, tukeywin, ...

%FT
SP = SP.FT(SP.data_stitched.t_stitched, permute(SP.data_stitched.data_R,[3 1 2]).*repmat(SP.window2.',[1 50 50])); % wavenumbers are in cm^-1, Raman spectrum is arbitrary units

SP = SP.make_raman_spectrum();



%% TEST 

            Total_delay=1/SP.Clock_Freq*SP.N_t*SP.DazzlerTimeConversion;
            time_axis=0:Total_delay/(SP.N_t-1):Total_delay;

            figure;
            subplot(2,3,1);
        
            hold on,
            plot(SP.data_stitched.t_stitched(:)*1E12,squeeze(mean(SP.data_stitched.data_R,[1 2])).*SP.window2.'); hold off
            xlabel('delay [ps]','fontsize',14);
            ylabel('Raman oscillations','fontsize',14);
            title('Integrated temporal response','fontsize',14);
            text(-0.1,1.1,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
            set(gca,'ytick',[])      

            subplot(2,3,2);
            plot(SP.wn,SP.ramanSpectrum);
            xlabel('Wavenumbers [cm^{-1}]','fontsize',14);
            ylabel('Raman Spectrum','fontsize',14);
            xlim([0,160])
            title('Integrated Raman Spectrum','fontsize',14)
            text(-0.1,1.1,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
            set(gca,'ytick',[]);

            subplot(2,3,3);
            imagesc(squeeze(mean(SP.data_stitched.data_R(:,:,:),3)));
            title('DC PD signal','fontsize',14); 
            xlabel('pixels');
            ylabel('pixels');
            colorbar;
            title('Integrated Raman Spectrum','fontsize',14)
            text(-0.1, 1.1,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)

            subplot(2,3,4);
            imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(1),:,:))));
            xlabel('pixels');
            ylabel('pixels');
            colorbar;
            title(sprintf('Image at %.1f cm^{-1}',SP.wn(SP.pixels_plot(1))));
            colormap('hot')
            
            subplot(2,3,5);
            imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(2),:,:))));
            xlabel('pixels');
            ylabel('pixels');
            colorbar;
            title(sprintf('Image at %.1f cm^{-1}',SP.wn(SP.pixels_plot(2))));
            
            subplot(2,3,6);
            imagesc(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(3),:,:))));
            xlabel('pixels');
            ylabel('pixels');
            colorbar;
            title(sprintf('Image at %.1f cm^{-1}',SP.wn(SP.pixels_plot(3))));

            figure(101),
            imgrgb = imagesc2rgb(squeeze(abs(SP.hyperspectralRamanImageComplex(SP.pixels_plot(1),:,:))));
            mask = roipoly(imgrgb);
            SP = SP.make_raman_spectrum_with_mask(mask);

            figure (102),
            plot(SP.wn,SP.ramanSpectrum);
            xlabel('Wavenumbers [cm^{-1}]','fontsize',14);
            ylabel('Raman Spectrum','fontsize',14);
            xlim([0,160])
            title('Integrated Raman Spectrum','fontsize',14)
            text(-0.1,1.1,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14)
            set(gca,'ytick',[]);
            % figure(101),clf,
            % for tt=1:5
            %     hold on,
            % 
            %     Raman=mean(SP.data_processed(tt).data_R,[2 3]);
            %     %T=mean(SP.data_processed(tt).data_T,[2 3]);
            %     plot(time_axis+SP.delays(tt)*1e-12,Raman);%./(1+T./max(T(:))))
            %    %hold on,plot(time_axis+SP.delays(tt)*1e-12,T/25,'bo')
            % 
            % end

% 
% 
% t_ini=repmat(time_axis,1,size(SP.delays,2));
% t_ini=t_ini(:);
% t_ini=t_ini(:).';
% for tt=1:length(SP.delays)-1
% indice=find(time_axis<(SP.delays(tt+1)-SP.delays(tt))*1e-12);
% t_delay_parts(tt,:)=SP.t0+SP.delays(tt)*1e-12+time_axis(indice);
% data_parts(tt,:)=SP.data_processed(tt).data_R(indice,:,:);
% end
% t_delay_parts(tt+1,:)=SP.t0+SP.delays(tt+1)*1e-12+time_axis(indice);
% data_parts(tt+1,:)=SP.data_processed(tt+1).data_R(indice,:,:);