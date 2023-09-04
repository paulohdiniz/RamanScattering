%% RAMAN IMPULSIVE EXCITATION :
clear variables
close all
clc
%% Load Azurite spectrum for example 
opts = delimitedTextImportOptions("NumVariables", 2);
% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ";";
% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2"];
opts.VariableTypes = ["double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
DefaultDataset2 = readtable("C:\Users\phdin\Desktop\PauloDiniz\CBZ_form_II.csv", opts);
A = table2array(DefaultDataset2);
A(:,2)=(A(:,2)-min(A(:,2)))/max(A(:,2));
clear opts
clear DefaultDataset2
%% LASER SPECTRUM
delta_lambda=0.2;
lambda_0=800;
lambda=(lambda_0-500):delta_lambda:(lambda_0+500);
BW=70; %nm
sigma=BW/2.3548;
spectrum=exp(-(lambda-lambda_0).^2./(2*sigma.^2));
figure,plot(lambda,(spectrum))
title('Spectrum')
% spectrum(spectrum<0.5)=0;
%% Raman spectrum
Raman_lines2=[20 45 70 85 95 120 170 190 225 500 1000]; % cm-1
Raman_lines=[5:50:8100];
for ii=1:length(lambda)
    for rr=1:length(Raman_lines)
        lambda_raman(ii,rr)=1/(1/lambda(ii)-Raman_lines(rr)/1e7);
    end
end
% create the sets of diracs :
for ii=1:length(lambda)
    for rr=1:length(Raman_lines)
        ind_raman=find(lambda-lambda_raman(ii,rr)>0,1);
        diracs(ii,rr,:)=[ii ind_raman];
    end
end
%% Plot on the spectrum
figure(1),clf,
plot(lambda,spectrum,'b','LineWidth',2)
rr=105;
ratio=0.5;
for ii=1:30:length(lambda)
    hold on,
    h2=plot([1 1]*lambda(diracs(ii,rr,1)),ratio*[0 1],'r','LineWidth',2);
    h3=plot([1 1]*lambda(diracs(ii,rr,2)),ratio*[0 1],'r','LineWidth',2);
    drawnow
    pause(0.1)
    delete(h2)
    delete(h3)
end
%
for ii=1:length(lambda)
    for rr=1:length(Raman_lines)
        excitation(ii,rr)=spectrum(diracs(ii,rr,1))*spectrum(diracs(ii,rr,2));
    end
end
%%
 Intensity=sum(excitation,1);
Intensity=Intensity/max(Intensity);
 figure(101),clf,%plot(Raman_lines,(Intensity),'linewidth',2)
for rr=1:length(Raman_lines2)
    ind_exc=find(Raman_lines-Raman_lines2(rr)>0,1);
    if Intensity(ind_exc)>0.06
        col='g';
        int=Intensity(ind_exc);
    else
        col='r';
        int=0.2;
    end
    hold on,plot([1 1]*Raman_lines2(rr),[0 int],col,'linewidth',2)
    hold on,plot([1 1]*Raman_lines2(rr),[0 int],col,'linewidth',2)
end
%% Randy's paper formula : 
tau_p=170*1e-15;
omega=linspace(0,3e15,1000);
D_omega=exp(-(tau_p*omega/2).^2/2);
figure(101),clf,hold on,
plot(omega/2/pi/3e8*1e-2,D_omega,'color',[00 208 000]/255,'linewidth',2)
hold on,plot(A(:,1),(A(:,2))/1.5,'b','linewidth',2)
visible_spectra=interp1(omega/2/pi/3e8*1e-2,D_omega,A(:,1));
hold on,plot(A(:,1),(A(:,2))/1.5.*visible_spectra,'r','linewidth',2)
xlim([0 3000])
set(gcf,'color','w')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',20)
xlabel('$\Omega$ (cm$^{-1}$) ','interpreter','latex','Fontsize',24)
ylabel('Intensity (u.a) ','interpreter','latex','Fontsize',24)