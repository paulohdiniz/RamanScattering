function [F_pos,Y_pos,Fs,L] = my_FT2_no_plot_Complex_nfft(time_vec,Data_vec,FTAxis,InverseOpt,power_of_2)
% Nicht plotten function
%%  Calculating the fourier transform, without fftshift
dt = diff(time_vec);
dt = dt(1:end-1);
dt = abs(dt);
dt = mean(dt);
Fs = 1/dt;
L = time_vec(end)-time_vec(1);
if (exist('FTAxis','var'))
     NFFT = 2^(nextpow2(size(Data_vec,FTAxis))+power_of_2);
    if InverseOpt
%         ReducedMeanSignal = bsxfun(@minus,Data_vec,mean(Data_vec,FTAxis));
        Y = (fft(Data_vec,NFFT,FTAxis));
        F = (((0:1/NFFT:1-1/NFFT)*Fs).');
        F_pos = fftshift(F-Fs/2);
        Y_pos = Y;
    else
        ReducedMeanSignal = bsxfun(@minus,Data_vec,mean(Data_vec,FTAxis));
        Y = dt*(fft(ReducedMeanSignal,NFFT,FTAxis));
        F = (((0:1/NFFT:1-1/NFFT)*Fs).');
        F_pos = fftshift(F-Fs/2);% By subtracting Fs/2, we first have all negative numbers, go through zero and then all positive numbers. Doing fftshift we move half of the coefficients from the end (this in our case are the positive frequencies) to the start of the array thereby giving us an array going from positive minimum to positive maximum and then have a negative maximum going up to negative minimum.
%         F_pos = fftshift(F-Fs/2);

        Y_pos = Y;
    end
%     B = shiftdim(Y,FTAxis-1);
%     C = B(1:floor(NFFT/2),:);
%     ArrayShape = size(B);
%     ArrayShape(1) = size(C,1);
%     Y_pos_shifted = reshape(C,ArrayShape);
%     Y_pos = shiftdim(Y_pos_shifted,length(size(Y_pos_shifted)-FTAxis));
else
    if InverseOpt
        NFFT = length(Data_vec);
        Y = (ifft(Data_vec,NFFT));
        F = (((0:1/NFFT:1-1/NFFT)*Fs).');

        F_pos = fftshift(F-Fs/2);
        Y_pos = Y;
    else
        NFFT = length(Data_vec);
        Y = (fft(Data_vec-mean(Data_vec),NFFT));
        F = (((0:1/NFFT:1-1/NFFT)*Fs).');
        F_pos = fftshift(F-Fs/2);
        Y_pos = Y;
    end
%     Y_pos = Y(1:floor(NFFT/2),:,:);

    
end


% figure;
% plot(F(1:NFFT/2),(Y(1:NFFT/2)));
% title('fourier transform');
% xlabel('Frequency (Hz)');
% ylabel('Fourier transform');
% grid;