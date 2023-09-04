close all
clc
clear variables
% Go where the image is
cd('C:\Users\Doctorant Sisira\Documents\')
%% Simulation of data from Scan image with Sid's experiment
% Generation of a damped Sinusoid. 
t=0:1e-6:1e-1; %time axis
f0=1000; %Hz
damping=f0/10;
tau0=5e-3;
ringing=sin(2*pi*f0*t+0.4*pi/2).*exp(-t*damping).*heaviside(t-tau0);
tau=1e-3
  CC=[0 10*diff(exp(-((t-tau0)/tau).^2))/1e-4/80];
figure,plot(t,ringing)
hold on,plot(t,CC)
hold on,plot(t,CC+0.2*ringing,'k')
signal=CC+0.2*ringing+0.2*rand(size(CC));
%%
n_breaks=3;
t_breaks=(1:(n_breaks))*max(t)/(n_breaks+1);
time_error=0.2*t_breaks(1)*rand(n_breaks,1);
t_stitch=t_breaks+time_error.';
for ii=1:n_breaks
    I_breaks(ii)=find(abs(t-t_breaks(ii))==min(abs(t-t_breaks(ii))));
    I_stitch(ii)=find(abs(t-t_stitch(ii))==min(abs(t-t_stitch(ii))));
%  figure(101),hold on,plot(abs(t-t_breaks(ii)),'o')
end
window_size=(I_breaks(2)-I_breaks(1));
 stitched_signal=[signal(1:I_breaks(1))];
 for ii=2:n_breaks+1
     if ii~=(n_breaks+1)
     stitched_signal=[stitched_signal signal(I_stitch(ii-1):I_stitch(ii-1)+window_size)];
     else
         stitched_signal=[stitched_signal signal(I_stitch(ii-1):end)];
         stitched_signal=[stitched_signal zeros(1,length(signal)-length(stitched_signal))];
     end
 end
 %%
figure,plot(t,signal,'--b')
hold on,plot(t,stitched_signal,'ro')
figure,plot(abs(fft(signal(8000:end))))
hold on,plot(abs(fft(stitched_signal(8000:end))))
xlim([0 200])
