
temp=SP.ramanSpectrum;
temp=temp./max(temp(:));

wn_temp=SP.wn(SP.wn<250);
temp=(temp(SP.wn<250));

wn_temp2=wn_temp(wn_temp>4);
temp=(temp(wn_temp>4));

[y,w,Q]=findpeaks(temp,wn_temp2,'MinPeakHeight',0.05,'MinPeakDistance',2,'MinPeakProminence',0.02)

figure,plot(wn_temp2,temp,'b')
hold on,plot(w,y,'ro')