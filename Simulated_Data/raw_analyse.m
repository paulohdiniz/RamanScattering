Total_delay=1/SP.Clock_Freq*SP.N_t*SP.DazzlerTimeConversion;
time_axis=(0:Total_delay/(SP.N_t-1):Total_delay)*1e12;

n_xp=1
dataR1 = squeeze(sum(sum(SP.data_raw(n_xp).data_R, 2), 3));
dataT1 = squeeze(sum(sum(SP.data_raw(n_xp).data_T, 2), 3));


dataR2 = squeeze(sum(sum(SP.data_processed(n_xp+1).data_R, 2), 3));
dataT2 = squeeze(sum(sum(SP.data_processed(n_xp+1).data_T, 2), 3));


Line_zidth=2
dataR1 = dataR1./max(dataR1(time_axis>2));
dataT1 = dataT1./max(dataT1);
range1=time_axis<3.5;
range2=time_axis>3.5;
figure(1),clf, plot(time_axis(range1), dataR1(range1),'b','Linewidth',Line_zidth), 
hold on,plot(time_axis(range2),dataR1(range2),'b--','Linewidth',Line_zidth)
hold on,plot(time_axis, dataT1,'r','Linewidth',Line_zidth)


dataR2 = dataR2./max(dataR2);
dataT2 = dataT2./max(dataT2);
 range3=time_axis<0.5;
 range4=time_axis>0.5;
figure(1),hold on, plot(time_axis(range4)+3, dataR2(range4),'color',[0 0.6 0],'Linewidth',Line_zidth),
hold on,plot(time_axis(range3)+3, dataR2(range3),'color',[0 0.6 0],'Linewidth',Line_zidth,'LineStyle','--'),
hold on,plot(time_axis+3, dataT2,'r','Linewidth',Line_zidth)
ylim([-1.5 1.5])

xlim([0 6])
set(gcf,'Color','w')
set(gca,"TickLabelInterpreter",'latex')
set(gca,"FontSize",24)
xlabel('Temps (ps)','Interpreter','latex','FontSize',24)
ylabel('Amplitude (u a)','Interpreter','latex','FontSize',24)

% 
% 
% sums = squeeze(sum(sum(data, 2), 3));
% sumsT = squeeze(sum(sum(dataT1, 2), 3));
% plot(sums, 'LineWidth', 0.5), hold on,
% plot(sumsT, 'LineWidth', 0.5), hold on,