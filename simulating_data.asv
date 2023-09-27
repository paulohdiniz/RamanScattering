close all
clc
clear variables
%% Simulation of data from sid's experiment to check the processing SP : 
% Global Parameters
c=3e8;
% parameters of Dazzler scans :
t_dazz=0:1.5000e-14:4.5e-12;
N_windows=5;
window_step=3e-12;
% Generation of the damped Sinusoids. 

w0=2*2*pi*[21 75 111]*1e2*c;
A0=[1 1 1];
Q_factor=[25 15 15]*5;
damping=w0./Q_factor;
% tau = [1 2 3];

tau0=1.7e-12;
% tau = [1 2 3]
for jj=0:N_windows-1
    time_temp=t_dazz+jj*window_step;
     [t_temp]=simulated_transmission(t_dazz,0.3e-12,40);
    for ii=1:length(w0)        
        ringing(ii,jj+1,:)=A0(ii)*sin(w0(ii)*time_temp+0.4*pi/2).*exp(-time_temp*damping(ii)).*heaviside(time_temp-tau0); %#ok<SAGROW>
         data_T(ii,jj+1,:)=t_temp; %#ok<SAGROW>
    end
      
end

%  [transmission]=simulated_transmission(t_dazz,0.3e-12,40);
%
tau0_cc=1.7e-12;

tau=0.2e-12;
  CC=[0 50*diff(exp(-((t_dazz-tau0_cc)/tau).^2))];
  figure,
for ii=1:N_windows
    hold on,
  plot((t_dazz+(ii-1)*3e-12)*1e12,squeeze(ringing(2,ii,:)))
  hold on,plot((t_dazz+(ii-1)*3e-12)*1e12,CC.*(ii==1))
  hold on,plot((t_dazz+(ii-1)*3e-12)*1e12,squeeze(sum(ringing(:,ii,:),1)).*squeeze(data_T(1,ii,:)),'k')
  hold on,plot((t_dazz+(ii-1)*3e-12)*1e12,squeeze(data_T(1,ii,:)),'r')
end


signal=ringing;
signal(1,1,:)=squeeze(signal(1,1,:)).'+CC;

% Noising ? 
 figure, plot(rand(301,1).*sqrt(abs(squeeze(sum(signal(1,1,:),1))*1000)).*sign(squeeze(sum(signal(1,1,:),1))))

%% Make the image
Nx=64; %size of the image
Ny=64;
gauss_filter_size=10; % size of the gaussian filter
img_temp=rand(Nx,Ny); % full random
img_temp(img_temp<.5)=0; % we threshold it a bit
img_temp(img_temp>.8)=1;% we threshold it a bit

img=medfilt2(img_temp,[gauss_filter_size,gauss_filter_size]); % First we gaussian filter it to get an image
img2=medfilt2(img,[gauss_filter_size,gauss_filter_size]); % We do it again to smooth the image
figure,imagesc((img2./max(img2(:))).^2)

%% Make the "scan image" image

img2_temp=repmat(img2./max(img2(:)),[1 1 length(t_dazz)]);

for ii=1:N_windows
            signal_temp=(sum(signal(:,ii,:),1));
    data_synth(ii).dataR=repmat(signal_temp,[Nx Ny 1]).*img2_temp;
%     data_synth(ii).dataT=repmat()
end

% THIS IS THE SYNTHETIC DATA

%%
signal=CC+0.2*ringing+0.2*rand(size(CC));
% %%
% n_breaks=3;
% t_breaks=(1:(n_breaks))*max(t)/(n_breaks+1);
% 
% time_error=0.2*t_breaks(1)*rand(n_breaks,1);
% 
% t_stitch=t_breaks+time_error.';
% 
% for ii=1:n_breaks
%     I_breaks(ii)=find(abs(t-t_breaks(ii))==min(abs(t-t_breaks(ii))));
%     I_stitch(ii)=find(abs(t-t_stitch(ii))==min(abs(t-t_stitch(ii))));
% %  figure(101),hold on,plot(abs(t-t_breaks(ii)),'o')
% end
% window_size=(I_breaks(2)-I_breaks(1));
%  stitched_signal=[signal(1:I_breaks(1))];
% 
%  for ii=2:n_breaks+1
% 
%      if ii~=(n_breaks+1)
%      stitched_signal=[stitched_signal signal(I_stitch(ii-1):I_stitch(ii-1)+window_size)];
%      else
%          stitched_signal=[stitched_signal signal(I_stitch(ii-1):end)];
%          stitched_signal=[stitched_signal zeros(1,length(signal)-length(stitched_signal))];
%      end
%  end
%  %%
% figure,plot(t,signal,'--b')
% hold on,plot(t,stitched_signal,'ro')
% 
% 
% figure,plot(abs(fft(signal(8000:end))))
% hold on,plot(abs(fft(stitched_signal(8000:end))))
% xlim([0 200])




