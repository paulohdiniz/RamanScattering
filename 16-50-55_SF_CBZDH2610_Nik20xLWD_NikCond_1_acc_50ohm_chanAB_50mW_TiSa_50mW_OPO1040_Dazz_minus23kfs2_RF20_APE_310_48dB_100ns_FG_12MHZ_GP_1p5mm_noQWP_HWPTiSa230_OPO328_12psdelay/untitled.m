 A=data_3d_parts_arranged;

  A=squeeze(sum(sum(A(:,:,1:50,:),3),4));
 A(1,1:90)=0;
figure(1),
 for ii=1:5
     hold on,plot(t_delay_parts_arranged(ii,:),A(ii,:)-mean(A(ii,:),2))
 end
 A=data_3d_parts_arranged;

   A=squeeze(sum(sum(A(:,:,50:end,:),3),4));

figure(1),
 for ii=1:5
     hold on,plot(t_delay_parts_arranged(ii,:),A(ii,:)-mean(A(ii,:),2))
 end