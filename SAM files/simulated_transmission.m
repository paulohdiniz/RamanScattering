function [transmission]=simulated_transmission (t_dazz,tau,dead_zone)
t_conv=[-t_dazz(end:-1:2) t_dazz];
temp=zeros(size(t_conv));
temp(length(t_dazz)+dead_zone:end-dead_zone)=1;

gauss=exp(-((t_conv)/tau).^2);
conv_temp=conv(temp,gauss);
% figure,plot(t_conv,gauss)
% hold on,plot(t_conv,temp)
% hold on,plot(t_conv,conv_temp(length(gauss):end))

t_temp=conv_temp(length(gauss):end);
t_temp=t_temp(1:length(t_dazz))./max(t_temp((1:length(t_dazz))));
noise=rand(1,6)/3;

noise=(interp1(linspace(dead_zone,length(t_dazz)-dead_zone,length(noise)),noise,dead_zone:length(t_dazz)-dead_zone,'cubic')).';
noise2=smooth([zeros(dead_zone,1).' noise.' zeros(dead_zone-1,1).'].',50);
figure,plot(smooth(t_temp+noise2.',20))
transmission=smooth(t_temp+noise2.',20);
end