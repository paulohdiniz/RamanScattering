close all;
clear all;
%%

DazzlerTimeConversion= 161E-15/1E-6;
DelayToStartFourier = 1.485e-12;%1.54 for CBZ-DH 1.339E-12
t0 = 0;
ReturnOpt = 0;
ZeroPadOpt = 0;
stitch = 1;
DCopt = 1;
DelayToEndFourier =15E-12;
SHGOpt =0;
%wns_plot=[21 75 111];% frequency in cm-1 %See line 232

%% When the measurement was saved using lab view, import the files first 

%In this case each measurement is a folder and each folder  contains files
%'0000.txt','0000.ascii' etc
pathname = uigetdir('/Users/paulohd/Desktop/PauloDiniz');
folders=dir(pathname);

folder_range=[15:19];%change this according to file used 
for ii=folder_range
cd(fullfile(pathname,folders(ii).name))
 data_filename{ii-folder_range(1)+1}=importdata('0000_.ASC');
 parameters_filename{ii-folder_range(1)+1}=importdata('0000_.txt');
% parametersXML_filename{ii-folder_range(1)+1}=importdata('0000_.XML');
% SHG_filename{ii-folder_range(1)+1} = importdata('0000__SHG.ASC');
end

%% Getting required parameters from diffrent files
for ii =folder_range
    cd(fullfile(pathname,folders(ii).name));
    folder=fullfile(pathname,folders(ii).name);
    dom=xmlread('0000_.XML');
    Find = dom.getElementsByTagName('Val');
    Val_ClockFreqValPix=Find.item(75);
    ClockFreqValPix=str2double(char(getNodeValue(Val_ClockFreqValPix.item(0))));%at first getNodeValue gives 'java.lang.String' typem convert to char and then double
    Val_DelayValPix=Find.item(75);
    Val_ScanRangeValPix=Find.item(103);
    ScanRangeValPix=str2double(char(getNodeValue(Val_ScanRangeValPix.item(0))));
  
    FindDelayValPix = (strfind(folder,'psdelay'));%in Sids code this value comes to 180 for data sid folder. Why I dont know. So for now I add a +1 to adjust for this and see if it works
     %should rewrite this part for the 2p5 delays

     
    if (ii-folder_range(1)+1) <=4 % I have written it this way because until 9ps delay we just haveto go one step to the left of the pixel containing 'psdelay'. The moment we have a two figure delay like 12ps then we have to include two pixels to the left.
        StitchDelay(ii-folder_range(1)+1) = str2double(folder(strfind(folder,'psdelay')-1))*1E-12;% here folder = 19-24-47_blackhair 100accum_20mW_TiSa_60mW_OPO_15cm_10cmlenses_48dBgain_240degphase_aomdrive_1.6V_optimizedcollectionobj_0psdelay_spot1_300nsTC_channelA+b_lengthanalysis. We look for the address of the beginng of the string 'psdelay'. strfind(folder,'psdelay')-1) goes one step back to see whether it is 0 or 3 or 6 ps. folder(strfind(folder,'psdelay')-1) gives us the string at the address strfind(folder,'psdelay')-1.  Then this string is converted to a numebr using str2double.
    elseif (ii-folder_range(1)+1)>4 
        StitchDelay(ii-folder_range(1)+1) = str2double(folder(FindDelayValPix-2:FindDelayValPix-1))*1E-12;
        %StitchDelay(i) = str2double(folder(strfind(folder,'psdelay')-1))*1E-12;
    end
    NumOfDelays =parameters_filename{ii-folder_range(1)+1}.data(11);
    NumOfPtsX = parameters_filename{ii-folder_range(1)+1}.data(1);
    NumOfPtsY = parameters_filename{ii-folder_range(1)+1}.data(2);

    
    if DCopt==1
        NumOfPtsX = NumOfPtsX*2;
    end
    
    c = 3e8; %[m/sec]
    TotDelayScan = 1/ClockFreqValPix * NumOfDelays * DazzlerTimeConversion; %[s]
    %TotDelayScan=4E-12;
    
    if(stitch)
        t_delay_parts(ii-folder_range(1)+1,:) = linspace(0,TotDelayScan,NumOfDelays) - t0 + StitchDelay(ii-folder_range(1)+1); %[s]
        x_delay = t_delay_parts * c * 100; %[cm]
        
        data_3d_reversed_time = reshape(data_filename{ii-folder_range(1)+1},[NumOfDelays,NumOfPtsX,NumOfPtsY]);
        data_3d_parts(ii-folder_range(1)+1,:,:,:) = flip(data_3d_reversed_time,1);% It is only here the flipping is happening
    else
        t_delay_parts = linspace(0,TotDelayScan,NumOfDelays) - t0;
        x_delay = t_delay_parts * c * 100;
        data_3d_reversed_time = reshape(data_filename{ii-folder_range(1)+1},[NumOfDelays,NumOfPtsX,NumOfPtsY]);
        data_3d = flip(data_3d_reversed_time,1);
        data_3d_ROI=data_3d;
        t_delay = t_delay_parts;
    end
    
end

   %% Arrange the stitched parts
    if (stitch)
        if(~DCopt)
            
            [StitchDelayArranged,StitchDelayArrangedIdx] = sort(StitchDelay);
            data_3d_parts_arranged = data_3d_parts(StitchDelayArrangedIdx,:,:,:);
            t_delay_parts_arranged = t_delay_parts(StitchDelayArrangedIdx,:);
            % The reason why he calls it t_delay_with_jumps is because at the
            % points where he is joining the first array (0ps delay) and the
            % fragmented (by cutting out the points that are already present in the
            % first array) second array(3ps delay) the successive time points are
            % not equally spaced.
            t_delay_with_jumps = squeeze(t_delay_parts_arranged(1,:));% the squeeze command kills the index(New dimension to account for the delayed measurements. i=1 for 0ps, 2 for 3 ps and so on) that we introduced to tag the element in sorted array(stich delay).
            data_3d_with_jumps = squeeze(data_3d_parts_arranged(1,:,:,:));
            for i=1:length(StitchDelayArrangedIdx(1:end-1))
                temp_t_vec = t_delay_parts_arranged(i+1,:);
                remove_t_indxs = temp_t_vec<=max(max(t_delay_with_jumps));% why TWO max commands? This step is to eliminate all entries in array 2 of time vector which are smaller than the max element of array 1. question is we compare with t_delay_with_jumps. It makes sense to compare an array from 3 to 7 ps with 4.5 ps. Why should we compare an array from 9 to 13 ps with 4.5 ps? Because according to this expression the max o of t_jump is always 4.25 ps
                temp_t_vec(remove_t_indxs) = [];% In this step we pass a logical array as an argument. We re killing cerain elements in the array that are repeated in the previous array.
                temp_data_vec = squeeze(data_3d_parts_arranged(i+1,:,:,:));% In this step we are trying to extract the 3D hyperspectral image acquired by setting a 3ps delay using the mechanical delay line. The 3D block is referred to by the value of the first index (i=2). We then squeeze the 4D array. This kills the dimension with singleton dimension which in this case is dimension which holds the index.
                temp_data_vec(remove_t_indxs,:,:) = [];%In the previous step we went from the 4D array to the 3D array.  Now we have passed to the first argument which is the delay axis an array which has the addresses of the delay values that were already accounted in the previous measurement.
                
                t_delay_with_jumps = cat(2,t_delay_with_jumps, squeeze(temp_t_vec));% Its a horizontal concatenation. What does squeeze here mean since we are applying it to a 1D array.
                data_3d_with_jumps = cat(1,data_3d_with_jumps, squeeze(temp_data_vec));% Its a vertical concatenation. the second squeeze here is for sanity check?
            end
        else
            [StitchDelayArranged,StitchDelayArrangedIdx] = sort(StitchDelay);
            data_3d_parts_arranged = data_3d_parts(StitchDelayArrangedIdx,:,:,:);
            t_delay_parts_arranged = t_delay_parts(StitchDelayArrangedIdx,:);
            
            
            t_delay_with_jumps = squeeze(t_delay_parts_arranged(1,:));
            data_3d_with_jumps = squeeze(data_3d_parts_arranged(1,:,:,:));%squeeze makes teh 4D matrix a 3D matrix. Here data_3d_with_jumps contains only the measurement corresponding to 0ps.
            Integrated_DCPD_parts_temp = squeeze(sum(sum(data_3d_with_jumps(:,NumOfPtsX/2+1:end,:),2),3));%Gives a 264*1 vector. to make it 1*264 we tranpose it later.
            for i=1:length(StitchDelayArrangedIdx(1:end-1))% here the index is i-1 because in the next step we take the i+1 th row of the array.
                temp_t_vec = t_delay_parts_arranged(i+1,:);% This assigns the time array from the successive time window. 
                overlap_t_indxs_curr = find(temp_t_vec<=max(max(t_delay_with_jumps)));
                overlap_t_indxs_prev = find(t_delay_with_jumps>=min(min(temp_t_vec)));
                
                t_ovrlp_curr = temp_t_vec(overlap_t_indxs_curr);
                
                PrevDCPDpxlsOvrlp = Integrated_DCPD_parts_temp(overlap_t_indxs_prev);%This is a column vector. Here we are comparing two integrated sums over the 50*50 pixels. The first sum considers 78 time steps at the end of the first measurement and the second sum considers 78 steps at the beginning of the second measurement.
                CurDCPDpxlsOvrlp = squeeze(sum(sum(data_3d_parts_arranged(i+1,overlap_t_indxs_curr,NumOfPtsX/2+1:end,:),3),4));%IS THERE A SPECIFIC REASON TO DO DEFINE IT AS A ROW VECTOR AND LATER TRANSPOSE IT? However, Because of the way we have used it here we have a 1*78 row vector. This is because we have not squeezed the data(1,:,:,:) 4D matrix before performing the sum. The first sum will give a 1*78*1*50 matrix. The second sum gives a 1*78. The final squeeze does nothing. Here the sum is over the 3rd and 4th indices. This is because we have not eaten the dummy index at the beginning(which labels the delay of acquisition. 1 for 0ps, 2 for 3 ps etc) with the squeeze function already.
                %CurDCPDpxlsOvrlp = squeeze(sum(sum(squeeze(data_3d_parts_arranged(1,overlap_t_indxs_curr,NumOfPtsX/2+1:end,:)),2),3));% This was justa trial made to see how the dimensions change when we sum across the 2 and 3rd index.
                %CurDCPDpxlsOvrlp = Integrated_DCPD_parts_temp(overlap_t_indxs_curr);% This works only when we stay in one time window and do not do any stitching in time.
                IntrsxnPix = find(abs(PrevDCPDpxlsOvrlp-CurDCPDpxlsOvrlp') == min(abs(PrevDCPDpxlsOvrlp-CurDCPDpxlsOvrlp')));%1) The prime on top of CurDCPDpxlsOvrlp' is to transpose the row vector to a column vector. This gets us to the intersection point of two DC curves in the successive delays.
                %IntrsxnPix = find(abs(PrevDCPDpxlsOvrlp-CurDCPDpxlsOvrlp) == min(abs(PrevDCPDpxlsOvrlp-CurDCPDpxlsOvrlp)));% Without prime.I have compensated for the prime by defining CurDCPDpxlsOvrlp through Integrated_DCPD_parts_temp. We find the difference between the two sums and try to get to the point at which the difference is minimum. This is the point at which the curves intersect.
                Intrsxn_T_curr = t_ovrlp_curr(IntrsxnPix);% The IntrsxnPix is found to be 52 for the HT69 sample.We pick one time value from the 78 time values that were common amongst the two windows (0ps amd 3ps). It is at this time value that the difference between the DC curves is minimum. 
                
                % we first try to kill all the elements beyond the
                % intersection point in the second time window
                temp_t_vec2 = temp_t_vec;
                temp_t_vec2(find(temp_t_vec<Intrsxn_T_curr)) = [];% This gives a 1*213 row vector. Out of 264 columns 51 are removed since we found the intersection pixel to be 52.
                temp_data_vec = squeeze(data_3d_parts_arranged(i+1,:,:,:));% to kill the singleton dimension i+1
                temp_data_vec(find(temp_t_vec<Intrsxn_T_curr),:,:) = [];%Kills the data recorded for all t's less than 3.8ps.
                
                
                data_3d_with_jumps = cat(1,data_3d_with_jumps(find(t_delay_with_jumps<Intrsxn_T_curr),:,:),temp_data_vec);% % Its a vertical concatenation. the second squeeze here is for sanity check?
                t_delay_with_jumps = cat(2,t_delay_with_jumps(find(t_delay_with_jumps<Intrsxn_T_curr)), squeeze(temp_t_vec2));% Here t_delay_with_jumps(find(t_delay_with_jumps<Intrsxn_T_curr)) gives a 1*237 vector which is concatenated with a 1*213 vector to give 1*450 vector.
                Integrated_DCPD_parts_temp = squeeze(sum(sum(data_3d_with_jumps(:,NumOfPtsX/2+1:end,:),2),3));%It is to be noted that  Integrated_DCPD_parts_temp is updated here since data_3d_with_jumps is concatenated for all the delays in the previous step. Here data_3d_with_jumps contains the measurements corresponding to both 0ps and 3 ps.
                
%                 figure;
%                 plot(t_delay_with_jumps,Integrated_DCPD_parts_temp);%To check how the concatenation happens
                
            end
                %t_delay_with_jumps=t_delay_with_jumps*1E12;% It is important to rescale t_delay_with_jumps because we will be interpolating t_delay with this. If we rescale here it is sufficient. The t_delay axis ill automatially be in ps. Earlier I thought this rescaling should be done in line

        end
        % Interpolate the parts to have equally spacing between all delay points (due to stitching)
        
        [interp_GridD,interp_GridX,interp_GridY] = ndgrid(t_delay_with_jumps,1:NumOfPtsX,1:NumOfPtsY);% We are setting up the grid. With t_delay_jumps on one axis, x points on one axis and y points on the other axis
        t_delay = linspace(min(t_delay_with_jumps),max(t_delay_with_jumps),length(t_delay_with_jumps));% Here we make the succesive points in the time vector equidistant so that we remove the jumps in the t_vec that were previous introduced when we concatenated time vectors corresponding to different delays. Here we make sure that the lengths of t_delay_with_jumps and t_delay are the same. Also we make sure that the beginning and the end of t_delay_with_jumps and t_delay are the same. This way we are more than good for interpolation.
        [interp_GridDFin,interp_GridXFin,interp_GridYFin] = ndgrid(t_delay,1:NumOfPtsX,1:NumOfPtsY);%Now we have a grid that has a equally spaced: time axis, x-axis and y-axis.
        
        data_3d = interpn(interp_GridD,interp_GridX,interp_GridY,data_3d_with_jumps,interp_GridDFin,interp_GridXFin,interp_GridYFin);%liner interpolation.  The data we have here is corrected for jumps. The time units are equally spaced. We have data corresponding to 0ps and 3ps stiched to each other.
        data_3d_ROI=data_3d;
    end  
    
   %% Remove offset and normalize by the DC value
    
    if(DCopt)
        
        %ModPD = data_3d(:,1:NumOfPtsX/2,:)  - squeeze(mean(mean(mean(data_3d_with_jumps(1:10,1:NumOfPtsX/2,:)))));%Data corresponding to the modulated part. Here we subtract the mean due to th first 10 points because we do not want to amplify and tiny fluctuation before zero delay. During these 10 steps the transmission of the dazzler is poor and when we normalize we might amplify any small non-zero value of the signal at these points. 5 steps = 80fs. So we consider the reading for the first 160 fs. By looking at the DC output plot I would say we can go twice this value. It is flat till 0.3 ps. The bottom line is that we dont want to amplify any small fluctuation before the zero delay. After the zero delay the transmission of the Dazzler is automatically larger and hence the normalization would not create any artefcat since it would only scale down the signal. Also it makes sense that when the transmission of the DAZZLER is low what we measure is the noise since the detector is measuring the 800nm which passes through the dazzler.
        ModPD = data_3d(:,1:NumOfPtsX/2,:)  - squeeze(mean(mean(mean(data_3d(1:10,1:NumOfPtsX/2,:)))));% I do it with data_3d and not data_3d_with_jumps because if I give stitch=1 and DCOpt =0 then I would not be generating data_3d_with jumps.Data corresponding to the modulated part. Here we subtract the mean due to th first 10 points because we do not want to amplify and tiny fluctuation before zero delay. During these 10 steps the transmission of the dazzler is poor and when we normalize we might amplify any small non-zero value of the signal at these points. 5 steps = 80fs. So we consider the reading for the first 160 fs. By looking at the DC output plot I would say we can go twice this value. It is flat till 0.3 ps. The bottom line is that we dont want to amplify any small fluctuation before the zero delay. After the zero delay the transmission of the Dazzler is automatically larger and hence the normalization would not create any artefcat since it would only scale down the signal. Also it makes sense that when the transmission of the DAZZLER is low what we measure is the noise since the detector is measuring the 800nm which passes through the dazzler.
        DCPD = data_3d(:,NumOfPtsX/2+1:end,:);% Data corresponding to the DC part
        
        DCPDTransmissionSum=squeeze(sum(sum(DCPD,2),3));
        DCPDTransmissionMean = squeeze(mean(mean(DCPD,2),3)); % This should give a 450*1 array. Each element of the array is a mean computed over all the spatial pixels(x and y) for a given delay.In the previous section we were performing a similar action where we computed the sum of the DC measurements which where not corrected for jumps.
        [DCPDTransmissionMean3D,~] = ndgrid(DCPDTransmissionMean,zeros(size(DCPD,2),size(DCPD,3)));% By doing this we construct a 450*2500 matrix. The 450*1 array is repeated 2500 times
        DCPDTransmissionMean3D = reshape(DCPDTransmissionMean3D,size(DCPD));% WE reshape the 450*2500 matrix into a 450*50*50 matrix. We do this so that we artifically build a matrix with the same dimensions as the data cube however with each pixel containing the mean transmission (as a function of the delay).
%         data_3d = ModPD./DCPDTransmissionMean3D;% Here we divide the time trace at each pixel (450*1) array by the integrated mean transmission (450*1 array).
        
        DC_output_normalized =  DCPDTransmissionMean/max(DCPDTransmissionMean);
        [DCPDTransmissionMean3D_norm,~] = ndgrid(DC_output_normalized,zeros(size(DCPD,2),size(DCPD,3)));% By doing this we construct a 450*2500 matrix. The 450*1 array is repeated 2500 times
        DCPDTransmissionMean3D_norm = reshape(DCPDTransmissionMean3D_norm,size(DCPD));% WE reshape the 450*2500 matrix into a 450*50*50 matrix. We do this so that we artifically build a matrix with the same dimensions as the data cube however with each pixel containing the mean transmission (as a function of the delay).
        data_3d_norm = ModPD./DCPDTransmissionMean3D_norm;% Here we divide the time trace at each pixel (450*1) array by the integrated mean transmission (450*1 array).
        data_3d = data_3d_norm; % i included this line on 23/8/19
        
    end
    
    
    
    %% Zero padding and Computing means and sums
    if ZeroPadOpt
        ZroPadFct = 2;
        integrated_time_trace = (squeeze(sum(sum(data_3d,2),3)));
        mean_time_trace = (squeeze(mean(mean(data_3d,2),3)));
        integrated_time_trace = padarray(integrated_time_trace,length(integrated_time_trace)*ZroPadFct,mean(integrated_time_trace(end-50:end)),'post');
        mean_time_trace = padarray(mean_time_trace,length(mean_time_trace)*ZroPadFct,mean(mean_time_trace(end-50:end)),'post');
        SummedImageOverDelay = squeeze(sum(abs(data_3d),1));%This is the red yellow orang pixelated image that we see on the bottom left corner. The squeeze here is effective since after the first sum we have a 1*50*50 array. squeeze eats the 1.
        MeanImageOverDelay = squeeze(mean(abs(data_3d),1));%This is the red yellow orang pixelated image that we see on the bottom left corner.
        dt = t_delay(2) - t_delay(1);
        t_delay_spect = min(t_delay):dt:(min(t_delay)+dt*(length(integrated_time_trace)-1));%this value can also be defined through mean_time_trace. We first pad the mean time trace with an array whose size is twice the size of the mean time trace array and the value of each entry of this padded array is the mean value of the last 50 elements of the mean time trace array.
    else
        integrated_time_trace = (squeeze(sum(sum(data_3d,2),3)));% Here the squeeze does nothing. The first sum returns a 450*1*50 array. The second sum returns a 450*1 array.
        mean_time_trace = (squeeze(mean(mean(data_3d,2),3)));% Here the squeeze does nothing. The first sum returns a 450*1*50 array. The second sum returns a 450*1 array.
        %integrated_time_trace_rescale=integrated_time_trace/max(integrated_time_trace);
        %integrated_time_trace = padarray(integrated_time_trace,5000,mean(integrated_time_trace(end-50:end)),'post');
        SummedImageOverDelay = squeeze(sum(abs(data_3d),1));%This is the red yellow orang pixelated image that we see on the bottom left corner. The squeeze here is effective since after the first sum we have a 1*50*50 array. squeeze eats the 1.
        MeanImageOverDelay = squeeze(mean(abs(data_3d),1));%This is the red yellow orang pixelated image that we see on the bottom left corner.
        t_delay_spect = t_delay;
    end
       
     %here
%% Find Fourier per pixel
    
    % PixelsForFourier = find(t_delay>=0);
    PixelsForFourier = find(t_delay>=DelayToStartFourier & t_delay<=DelayToEndFourier);
    PixelsForFourierSpect = find(t_delay_spect>=DelayToStartFourier & t_delay_spect<=DelayToEndFourier);

    x_delay = c*t_delay * 100;
    x_delaySpect = c*t_delay_spect * 100;
    
    
%     [wn,RamanSpectrum] = my_FT2_no_plot(x_delaySpect(PixelsForFourierSpect),integrated_time_trace(PixelsForFourierSpect)); % wavenumbers are in cm^-1, Raman spectrum is arbitrary units
%     [wn,HyperspectralRamanImage] = my_FT2_no_plot(x_delay(PixelsForFourier),data_3d(PixelsForFourier,:,:),1); % wavenumbers are in cm^-1, Raman spectrum is arbitrary units    
   
    [wn,HyperspectralRamanImageComplex] = my_FT2_no_plot_Complex_nfft(x_delay(PixelsForFourier), data_3d(PixelsForFourier,:,:),1,0,1); % wavenumbers are in cm^-1, Raman spectrum is arbitrary units
    wn = wn(1:floor(size(HyperspectralRamanImageComplex,1)/2));
    %     figure; plot(tempwn,abs(sum(sum(HyperspectralRamanImageComplex,2),3)));
    HyperspectralRamanImageComplex(abs(wn)>900,:,:) = 0;
%     [tempt,fourierFilteredSig] = my_FT2_no_plot_Complex(tempwn,HyperspectralRamanImageComplex,1,1);
   %  figure; plot((real(sum(sum(fourierFilteredSig,2),3))));

    PosFreqHyperspectralRamanImageComplex = HyperspectralRamanImageComplex(1:floor(size(HyperspectralRamanImageComplex,1)/2),:,:);
    HyperspectralRamanImage = abs(PosFreqHyperspectralRamanImageComplex);
  
   RamanSpectrum = abs(sum(sum(PosFreqHyperspectralRamanImageComplex,2),3));

   %Here it takes the values of the 3 biggest peaks and saves it in the wns_plot variable to plot later
   [yPeaks,xPeaks] = findpeaks(RamanSpectrum,wn,'SortStr','descend');
   wns_plot=[xPeaks(1) xPeaks(2) xPeaks(3)];
   

    %% To compute SNR
    
    % Defining the pixel at which we consider the noise. That is, We choose
    % the pixel at which we know for sure that there is no sample and find
    % the variance in this pixel.
   
    Noise_pixel_array = ModPD(:,9,5);%ModPD(:,9,10)
    Std_noise=std(Noise_pixel_array);
    
    %To see how the noise looks like
    %
    
    % Defining the pixel at which we consider the signal. That is, We choose
    % the pixel at which we know for sure that there is sample and find
    % the rms value (mean square since we are measuring the voltage) in this pixel.
    
    signal_pixel_array = ModPD(:,19,14);%(31,25)
    
    mean_signal = mean(signal_pixel_array);
    square_signal = (signal_pixel_array).^2;
    mean_square_signal=mean(square_signal);
%     
%     %To see how the signal looks like
%     figure;set(gcf,'Color','white');fig.Position=[0 50 1850 900];
%     
%     plot(t_delay,signal_pixel_array);
%     hold on;
%     plot(t_delay,Noise_pixel_array);
    
    % The GOLDEN NUMBER
    SNR_through_variance =  mean_square_signal/(Std_noise)^2;
    SNR_through_variance_dB = 10*log(SNR_through_variance);
    SNR_through_std =  (mean_signal)/(Std_noise);
   %% To plot stuff

TitleFontSize = 14;
AxisFont = 14;
SubplotLabelPos(1) = -0.1;
SubplotLabelPos(2) = 1.1;
ImageFoV=200;
ScaleBarLength = (50/ImageFoV)*2; %NumOfPxls/ImageFoV_um*NumOf_um
%%

for i=1:size(wns_plot,2)
    pixels_plot(i) = find(abs(wn-wns_plot(i))==min(abs(wn-wns_plot(i))));
end
% The next three lines dont affect the rest of the code
WnPix = pixels_plot(1);
MeanImage_wn1 = squeeze(mean(HyperspectralRamanImage(WnPix-1:WnPix+1,:,:),1));%
Max_int_pixel = find(MeanImage_wn1==max(max(MeanImage_wn1)));
Max_int_pixel = ceil(Max_int_pixel/50);
% WnPix = pixels_plot(2);
% MeanImage_wn2 = mean(HyperspectralRamanImage(WnPix-1:WnPix+1,:,:),1);%
% Max_int_pixel = max(max(MeanImage_wn2,1),2);
% 
% WnPix = pixels_plot(3);
% MeanImage_wn3 = mean(HyperspectralRamanImage(WnPix-1:WnPix+1,:,:),1);%
% Max_int_pixel = max(max(MeanImage_wn3,1),2);
ThrshVal = 4*std(RamanSpectrum(find(wn>160 & wn<250)));
PixelsAboveThreshold = find(RamanSpectrum>ThrshVal);
MeanRamanImage = mean(HyperspectralRamanImage(PixelsAboveThreshold,:,:),1);
% The MeanRamanImage has not been used elsewhere. 
% I feel it is non sensical to even have this quantity because 
% we are averaging over a range of wave numbers at which the signal is greater than 4std. 
% This neither gives the integrated time trace nor gives the Raman image.

%% To plot without scale bar

figure;
subplot(2,3,1);
plot(t_delay_spect*1E12,integrated_time_trace);hold on
plot(t_delay_spect(PixelsForFourierSpect(1))*1E12,integrated_time_trace(PixelsForFourierSpect(1)),'o'); hold off
%plot(t_delay_spect(PixelsForFourierSpect(1))*1E12,integrated_time_trace(PixelsForFourierSpect(1)),'o'); hold off
xlabel('delay [ps]','fontsize',AxisFont);
ylabel('Raman oscillations','fontsize',AxisFont);
title('Integrated temporal response','fontsize',TitleFontSize)
ylim(1.2*[min(min(integrated_time_trace(PixelsForFourierSpect))), max(max(integrated_time_trace(PixelsForFourierSpect)))]);
xlim([0 9]);
text(SubplotLabelPos(1),SubplotLabelPos(2),'a','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
set(gca,'ytick',[])


 subplot(2,3,2);
 plot(wn,RamanSpectrum);
 xlabel('Wavenumbers [cm^{-1}]','fontsize',AxisFont);
 ylabel('Raman Spectrum','fontsize',AxisFont);
 xlim([0,160])
 title('Integrated Raman Spectrum','fontsize',TitleFontSize)
 text(SubplotLabelPos(1),SubplotLabelPos(2),'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
 set(gca,'ytick',[]);
 
 
 subplot(2,3,3);
 if DCopt
     imagesc(squeeze(mean(DCPD,1)));
     title('DC PD signal','fontsize',TitleFontSize);
 else
     imagesc(SummedImageOverDelay);
     title('Integration over all delays','fontsize',TitleFontSize);
 end
 xlabel('pixels');
 ylabel('pixels');
 colorbar;
 title('Integrated Raman Spectrum','fontsize',TitleFontSize)
 text(SubplotLabelPos(1),SubplotLabelPos(2),'c','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
 
 WnPix = pixels_plot(1);
 subplot(2,3,4);
 imagesc(squeeze(HyperspectralRamanImage(WnPix,:,:)));
 xlabel('pixels');
 ylabel('pixels');
 colorbar;
 title(sprintf('Image at %.1f cm^{-1}',wn(WnPix)));
 colormap('hot')
 
 
 WnPix = pixels_plot(2);
 subplot(2,3,5);
 imagesc(squeeze(HyperspectralRamanImage(WnPix,:,:)));
 xlabel('pixels');
 ylabel('pixels');
 colorbar;
 title(sprintf('Image at %.1f cm^{-1}',wn(WnPix)));
 
 WnPix = pixels_plot(3);
 subplot(2,3,6);
 imagesc(squeeze(HyperspectralRamanImage(WnPix,:,:)));
 xlabel('pixels');
 ylabel('pixels');
 colorbar;
 title(sprintf('Image at %.1f cm^{-1}',wn(WnPix)));

%% Extraction from spontaneous raman

% if anthracene
% open('SpontaneousAntraceneRaman.fig');
% a = get(gca,'Children');
% xdata = get(a, 'XData');
% ydata = get(a, 'YData');
% zdata = get(a, 'ZData');
% 
% xdata3 = xdata(510:785);
% ydata3 = ydata(510:785);
% end
% 
% if CBZDH
% open('CBZDH.fig');
% a = get(gca,'Children');
% xdata = get(a, 'XData');
% ydata = get(a, 'YData');
% zdata = get(a, 'ZData');
% 
% xdata3 = xdata(510:785);
% ydata3 = ydata(510:785);
% end
       
%% To plot with scale bar
% This is the SNR in the particular pixel that I have defined earlier

fig=figure ;set(gcf,'Color','white');fig.Position=[0 50 1850 900];
plot(t_delay*1E12,signal_pixel_array);
hold on;
plot(t_delay*1E12,Noise_pixel_array);
xlabel('delay [ps]','fontsize',AxisFont);
ylabel('Intensity change','fontsize',AxisFont);
title(['Signal vs noise, SNR = ', num2str(SNR_through_variance)],'fontsize',TitleFontSize);
legend('Signal','Noise');
% xlim([0 16]);
ylim(1.2*[min(signal_pixel_array(PixelsForFourierSpect)) max(signal_pixel_array(PixelsForFourierSpect))]);

 fig=figure; set(gcf,'Color','white');fig.Position=[0 50 1850 900];
 subplot(2,3,1);
 plot(t_delay_spect*1E12,mean_time_trace);hold on
 plot(t_delay_spect(PixelsForFourierSpect(1))*1E12,mean_time_trace(PixelsForFourierSpect(1)),'o'); hold off% This line is just to have the marker
 xlabel('delay [ps]','fontsize',AxisFont);
 ylabel('Intensity change','fontsize',AxisFont);
 title('Integrated temporal response','fontsize',TitleFontSize)
 ylim(1.2*[min(min(mean_time_trace(PixelsForFourierSpect))), max(max(mean_time_trace(PixelsForFourierSpect)))]);
  xlim([0 15]);
 %ylim([0 5]);
 text(SubplotLabelPos(1),SubplotLabelPos(2),'a','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
 set(gca,'ytick',[])
%  legend('Noise','Signal')
 
 subplot(2,3,2);
 
%  yyaxis left;
%  plot(wn,RamanSpectrum);
%  xlabel('Wavenumbers [cm^{-1}]','fontsize',AxisFont);
%  ylabel('Amplitude','fontsize',AxisFont);
%  set(gca,'ytick',[]);
%  
%  %COMMENT THE LINES BELOW IF THERE IS NO SPONTANEOUS PLOT
%  hold on;
%  
%  yyaxis right;
%  plot(xdata3,ydata3);
%  xlabel('Wavenumbers [cm^{-1}]','fontsize',AxisFont);
%  %ylabel('Spontaneous Raman spectrum (a.u)','fontsize',AxisFont);
%  set(gca,'ytick',[]);
%  
 
%  xlim([0,200])
%  title('Integrated Raman Spectrum','fontsize',TitleFontSize)
%  text(SubplotLabelPos(1),SubplotLabelPos(2),'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
 
 
 SpntData = extract_data_from_2d_graph(SpontDataAdress);
        plot(wn,RamanSpectrum); hold on ;
        plot(SpntData(1,:),SpntData(2,:)/max(SpntData(2,:))*max(RamanSpectrum),'r-'); hold off
                %plot(wn,ones(size(wn))*ThrshVal,'r-'); hold off
        xlabel('Wavenumbers [cm^{-1}]','fontsize',AxisFont);
                ylabel('Raman Spectrum','fontsize',AxisFont);
        xlim([0,160])
        title('Integrated Raman Spectrum','fontsize',TitleFontSize)
        text(SubplotLabelPos(1),SubplotLabelPos(2),'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
%         set(gca,'ytick',[])
 
 
 subplot(2,3,3);
 if DCopt
     MeanImage = mean(DCPD,1);% Changing the range to 5:4+ScaleBarLength instead of 5:5+ScaleBarLength because we end up having one extra pixel. If scale bar length is 10 it means we wish to have 10 pixels in the scale bar. Going from 5:5+scalebarlength we have 5:15 which is 0 to 10 = 11 pixels. 
     MeanImage(:,46:47,5:4+ScaleBarLength) = max(max(MeanImage));%Its important to have this line before imagesc because we are actually writing on the color map of the actual image
     imagesc(squeeze(MeanImage));
     title('DC PD signal','fontsize',TitleFontSize);
 else
     imagesc(SummedImageOverDelay);
     title('Integration over all delays','fontsize',TitleFontSize);
 end
 colormap('hot');
 xlabel('pixels');
 ylabel('pixels');
 colorbar;
 %title('Integrated Raman Spectrum','fontsize',TitleFontSize)
 text(SubplotLabelPos(1),SubplotLabelPos(2),'c','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
 axis off
%  set(colorbar,'YTick',[])
 
 WnPix = pixels_plot(1);
 subplot(2,3,4);
  MeanImage = mean(HyperspectralRamanImage(WnPix-1:WnPix+1,:,:),1);%
 MeanImage(:,46:47,5:4+ScaleBarLength) = max(max(MeanImage));% Changing the range to 5:4+ScaleBarLength instead of 5:5+ScaleBarLength because we end up having one extra pixel. If scale bar length is 10 it means we wish to have 10 pixels in the scale bar. Going from 5:5+scalebarlength we have 5:15 which is 0 to 10 = 11 pixels. 
 imagesc(squeeze(MeanImage));
 xlabel('pixels','fontsize',AxisFont);
 ylabel('pixels','fontsize',AxisFont);
 colorbar;
 title(sprintf('Image at %.0f cm^{-1}',wn(WnPix)),'fontsize',TitleFontSize);
 text(SubplotLabelPos(1),SubplotLabelPos(2),'d','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
 colormap('hot')
 axis off
 axis square
%  set(colorbar,'YTick',[])
 
 WnPix = pixels_plot(2);
 subplot(2,3,5);
 MeanImage = mean(HyperspectralRamanImage(WnPix-1:WnPix+1,:,:),1);%
 MeanImage(:,46:47,5:4+ScaleBarLength) = max(max(MeanImage));% Changing the range to 5:4+ScaleBarLength instead of 5:5+ScaleBarLength because we end up having one extra pixel. If scale bar length is 10 it means we wish to have 10 pixels in the scale bar. Going from 5:5+scalebarlength we have 5:15 which is 0 to 10 = 11 pixels. 
 imagesc(squeeze(MeanImage));
 xlabel('pixels','fontsize',AxisFont);
 ylabel('pixels','fontsize',AxisFont);
 colorbar;
 title(sprintf('Image at %.0f cm^{-1}',wn(WnPix)),'fontsize',TitleFontSize);
 text(SubplotLabelPos(1),SubplotLabelPos(2),'e','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
 colormap('hot')
 axis off
 axis square
%  set(colorbar,'YTick',[])
 
 WnPix = pixels_plot(3);
 subplot(2,3,6);
 MeanImage = mean(HyperspectralRamanImage(WnPix-1:WnPix+1,:,:),1);%
 MeanImage(:,46:47,5:4+ScaleBarLength) = max(max(MeanImage));% Changing the range to 5:4+ScaleBarLength instead of 5:5+ScaleBarLength because we end up having one extra pixel. If scale bar length is 10 it means we wish to have 10 pixels in the scale bar. Going from 5:5+scalebarlength we have 5:15 which is 0 to 10 = 11 pixels. 
 if SHGOpt
     MeanImage = SHGdata;
     MeanImage(46:47,5:4+ScaleBarLength) = max(max(MeanImage));
 end
 imagesc(squeeze(MeanImage));
 xlabel('pixels','fontsize',AxisFont);
 ylabel('pixels','fontsize',AxisFont);
 title(sprintf('Image at %.0f cm^{-1}',wn(WnPix)),'fontsize',TitleFontSize);
 text(SubplotLabelPos(1),SubplotLabelPos(2),'f','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
 if SHGOpt
 title(sprintf('SHG image'),'fontsize',TitleFontSize);
 end
 colorbar;
 colormap('hot')
 axis off
 axis square
%  set(colorbar,'YTick',[])
   
 fig=figure; set(gcf,'Color','white');fig.Position=[0 50 1850 900];
 plot(wn,RamanSpectrum);
 %         plot(SpntData(1,:),SpntData(2,:)/max(SpntData(2,:))*max(RamanSpectrum),'r-'); hold off
 %plot(wn,ones(size(wn))*ThrshVal,'r-'); hold off
 xlabel('Wavenumbers [cm^{-1}]','fontsize',AxisFont);
 ylabel('Raman Spectrum','fontsize',AxisFont);
 xlim([0,160])
 title('Integrated Raman Spectrum','fontsize',TitleFontSize)
 text(SubplotLabelPos(1),SubplotLabelPos(2),'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',AxisFont)
%  set(gca,'ytick',[])