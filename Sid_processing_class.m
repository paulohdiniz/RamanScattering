close all
clear variables 
cd('/Users/paulohd/Desktop/PauloDiniz')
addpath(genpath('D:\Matlab'))
addpath('/Users/paulohd/Desktop/PauloDiniz')
%% CODE FOR USING THE CLASS THAT MAKES THE SAME THING AS SID'S CODE
% We generate the class
SP=Sid_Processing();
SP.DCopt=1;
% Choose the experiment, and load all the data from the correct folders
SP=SP.choose_folders_load_data(3);

tukey_window_param=0.1;
deadtime=74;
SP=SP.window_overlap(tukey_window_param,deadtime);
% SP=SP.Tnorm_and_center_data(1,0,deadtime);

  SP=SP.stitch_time_axis_T_with_interp('makima')

%% TEST 

     Total_delay=1/SP.Clock_Freq*SP.N_t*SP.DazzlerTimeConversion;
            time_axis=0:Total_delay/(SP.N_t-1):Total_delay;

            figure(101),clf,

            for tt=1:4
                hold on,
                Raman=mean(SP.data_processed(tt).data_R,[2 3]);
                T=mean(SP.data_processed(tt).data_T,[2 3]);
               plot(time_axis+SP.delays(tt)*1e-12,Raman);%./(1+T./max(T(:))))
               hold on,plot(time_axis+SP.delays(tt)*1e-12,T/25,'bo')
              
            end

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