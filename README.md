
## Overview

This repository contains MATLAB scripts, classes and functions for processing impulsive Raman diffusion data stimulated by the spectral decaling method. The processing involves various steps, from loading and normalizing data to generating Raman spectra and optimizing image quality. The code is designed to work with data files generated by LabView in the form of .ASC, .XML, and .TXT files.

These files are organized into five folders, sharing identical names except for the suffix indicating the spectral delay. As a result, our data processing involves extracting information from folders of the following structure: 

> &#128193; KL_CBZDH2610_Nik20xLWD_NikCond_1_acc_50ohm_chanAB_50mW_TiSa_50mW_OPO1040_Dazz_minus23kfs2_RF20_APE_310_48dB_100ns_FG_12MHZ_GP_1p5mm_noQWP_HWPTiSa230_OPO328_**0psdelay**

Summary of the lab files can be found in the *DetailsOfArchives* file. The first column indicates the date, the second the type of experiment, the third the number of experiments carried out, the fourth column indicates whether some data were changed or not and the fifth column a general summary of the quality of the images.

This repository contains two samples. The first is a sample of carbamazepine dihydrate (CBZDH) alone and the other is a sample containing CBZDH and anthracene. In our catalog they are samples 210617-xp4 and 210507-xp41 respectively.

## Requirements

Ensure you have MATLAB installed, preferably version **R2023a** or a more recent release, with the following toolboxes:

- *Computer Vision Toolbox, version 10.4*
- *Curve Fitting Toolbox, version 3.9*
- *Image Acquisition Toolbox, version 6.7.1*
- *Image Processing Toolbox, version 11.7*
- *RF Toolbox, version 4.5*
- *Signal Integrity Toolbox, version 13*
- *Signal Processing Toolbox, version 9.2*
- *Statistics and Machine Learning Toolbox, version 12.5*
- *Symbolic Math Toolbox, version 9.3*

## Getting Started

1. Clone this repository to your local machine.

2. Open MATLAB and navigate to the repository directory.

3. Set the path to the folder containing the five folders to be analyzed. Example:

    ```matlab
    dataPath = 'C:\Users\phdin\Desktop\Data_Paulo\210507';
    ```

4. Open the `Raman_processing_main.m` and put your path in the cd(path) and addpath(path) line. Choose the number of the experiment you want to work on, if you have 50 folders within your path corresponding to 10 experiments, and you want to choose experiment 7, change the xp_number to 7.

## Organization 

To process post-experiment data, we created the `Raman_Processing.m` class, which controls all the steps necessary for processing process data. The `Raman_processing_main.m` file is the script that executes what we want to do in this treatment.

To optimize the code, we use the `RmFast.m` class and the `RmFast_Main.m` script. The performance of code execution is essential so that we can work at the video rate. This code is a shorter code, but it uses the results of the analysis carried out in the previous data processing code. It is a code that tries to optimize the time but also the quality of the Raman image.

We'll start with the code for post-treatment, running the `Raman_processing_main.m` script will create an object called RP that will be updated when each method is executed.

## 1. Data Loading 

The first method is `choose_folders_load_data`. Firstly, it takes the name of the folders (like the example above) and gets from them the characteristics and parameters of the *Function Generator*, *Lockin*, *Objectives*, *LabView*, *Laser*, *Dazzler* and *Polarization*. From the .ASC file it will take the raw data, and other experiment parameters will be taken from the .TXT and .XML files.

## 2. Preprocessing
  ### 2.1 Eliminating Initial Impulse
  After that, the script executes the `window_overlap` method which will take *tukey_window_param* and *percent_FWHM* as parameters. This function serves to remove the large value of the signal caused by the initial impulse, therefore, it will only be removed from the first folder that contains data from 0ps to ~3ps. The *percent_FWHM* value is the percentage used of the full width at half maximum value of the square of the absolute value of the raw curve integral. We found that the 100% value leaves the value of the initial impulse insignificant in relation to the oscillation produced in the studied species, however, a window showing all these curves will open to the user to ratify or change this value. The *tukey_window_param value* indicates the percentage of the [Tukey window](https://fr.mathworks.com/help/signal/ref/tukeywin.html) used to smooth this initially made cut.
  ### 2.2 Normalization and centralization
  After that, the script executes the `Tnorm_and_center_data` method that normalizes and centers the data.
  ### 2.3 Stitching and interpolating
  The `stitch_time_axis_T_with_interp` method interpolates the data from the first folder with the initial impulse cut and the four other folders. For the issue of image quality, the interpolation done by the [*makima*](https://fr.mathworks.com/help/matlab/ref/makima.html) method is better done, but for performance it is the slowest, therefore it will be used in post-processing analysis but not in the performance code. Other methods can be used such as: [nearest](https://fr.mathworks.com/help/matlab/ref/interp1.html), [linear](https://fr.mathworks.com/help/matlab/ref/interp1.html), [pchip](https://fr.mathworks.com/help/matlab/ref/pchip.html) and [spline](https://fr.mathworks.com/help/matlab/ref/spline.html).
  ### 2.4 Windowing
  After interpolating the entire signal, we apply a window across its entire length. We can select the window with the `pick_fourier_window` method. MATLAB provides over [16 window options](https://fr.mathworks.com/help/signal/ug/windows.html), all of which are also accessible in our program. Users have the flexibility to select the specific ratio they wish to employ for this window.
  
## 3. Fourier Transform, Raman Spectrum and Hyperspectral Images

Following this step, the Fourier transform is executed using the `FT` method, which is performed using the [FFT](https://fr.mathworks.com/help/matlab/ref/fft.html) (Fast Fourier Transform) function. Subsequently, the Raman spectrum is generated using the `make_raman_spectrum` method.

## 4. Collecting image quality data

After the previous step, a data cube with the hyperspectral images is created. At this point we can analyze each image related to a certain wave number. Using the `calculated_images_scores_per_wn` method, we collect information about image quality methods, use the [*SSIM*](https://fr.mathworks.com/help/images/ref/ssim.html) (with the transmission image as reference) and the unreferenced methods [*piqe*](https://fr.mathworks.com/help/images/ref/piqe.html), [*niqe*](https://fr.mathworks.com/help/images/ref/niqe.html), and [*brisque*](https://fr.mathworks.com/help/images/ref/brisque.html), and save these values for later analyse.

## 5. Peak Detection


Subsequently, the `points_to_plot_by_frequency` method utilizes the [findpeaks](https://fr.mathworks.com/help/signal/ref/findpeaks.html) function to identify the peaks in the Raman plot by [prominence](https://fr.mathworks.com/help/signal/ug/prominence.html). This search is defined by the user who can stipulate the range sought and the minimum values between peaks, etc., the standard attributes are: 

```matlab
% Peak Detection Configurations
Max_wn_to_search_peaks = 250;      % Maximum wavenumber to search for peaks.
Min_wn_to_search_peaks = 5;        % Minimum wavenumber to search for peaks.
MinPeakHeight = 0.05;              % Minimum peak height (5% of the maximum).
MinPeakDistance = 10;              % Minimum distance between peaks in cm^-1.
MinPeakProminence = 0.02;          % Minimum peak prominence (2% of the maximum).

```

## 6. Inverse Fourier Transform

The last method is `get_signal_by_time_from_ifft` which, based on the peaks found, performs the [inverse fourier transform](https://fr.mathworks.com/help/matlab/ref/ifft.html) to find the signal in time and stores the oscillations of the peaks in the *signalIFFT* variable.

## FUNCTIONS
- **`plot_graphs`**: is the function that shows the graph of Raman oscillations over time (in ps) and frequency (in cm^-1) and alongside the three images corresponding to the peaks marked in the frequency spectrum. The time graph shows the raw signal (after interpolation), the window that will be used in it and the interpolated signal after using the window. The SSIM values associated with each image are relative to the transmission image.

- ![plot graphs image.](/Images/01.svg)

- **`plot_graphs_with_transmission`**: It has the same graphics but this time with the addition of the transmission image.

- ![plot graphs image.](/Images/02.svg)

- **`plot_graphs_with_roi`**: the same previous image opens, but this time a second image opens alongside with the transmission image, in this image the user can select two regions of interest to analyze the spectrum of the chosen region. First, the user presses the mouse with the left button, selecting points in the image and places the number of points he wants until the polygon is closed on the first selected button, at which point the user clicks twice to "send" the region for analysis and redoes a second polygon with the second area for analysis, after the second sending a third window will open showing the total spectrum and the spectrum of each selected region.

- ![roi1 gif.](/Images/roi1.gif)

- **`plot_similar_pixels_from_rois`**: this function shows the transmission image to the user who can select two areas of interest, exactly the same as the previous function, however, here the user can only take a part where he is sure there is substance A and then he selects only a part where it knows that substance B exists, so this function using MATLAB's correlation function will show the user the pixels where there is a similarity between their curves, this was tested and proved to be very efficient for our samples.

- ![roi2 gif.](/Images/roi2.gif)

- **`plot_graphs_with_ifft`**: next to the time and frequency graphs the user will have three curves corresponding to the inverse fourier transform for the selected peaks. (The time axis is not yet well adjusted).

- ![plot graphs image.](/Images/05.svg)

- **`plot_best_ssim_by_ratio_window`**: This function will give two graphs, one graph of the average of the ssims of the three images by the ratio_window and the other is the Raman spectrum according to the variation of the ratio_window, the darker the closer to 1. The function is very slow, and was made just for an interval of 0.1 of ratio_window, for more detailed curves and a higher density of points, it is necessary to decrease the step to a smaller value. However, this adjustment will require more time to generate the two graphs.

- ![plot graphs image.](/Images/06.svg)

- **`plot_spectrogram`**: A window will open showing the pspectrum, persistence spectrum and power spectrum graphs, the user must change the observation limit frequency..

- ![plot graphs image.](/Images/07.svg)

- **`plot_best_hyperspectral_images`**: in the previous spectrogram we can see that some curves end their signal before the end of the graph (~16 ps), so using the data until the end to build the image causes us to add unnecessary noise, so this function takes this ideal time where the curve enters the noise level and plots the spectrum image going only to that ideal time. If the image does not fall to the noise level, we use all the time to build the image, however, if it falls earlier we only use that time to optimize the image quality. The window formed by this function has in its first column the data used with the respective ideal time in ps, the second column indicates which line of the spectrum is being analyzed, the third column is the side view of this line collected from the spectrogram, the fourth is the image with as little added noise as possible.

- ![plot graphs image.](/Images/08.svg)

- **`plot_images_with_filters_by_psnr`**: This function takes the transmission image and looks for the three best filters (which give the highest psnr) and writes them into columns. So it takes the three ideal images that we saw in the previous function and uses these filters on them and uses the psnr in relation to before and after using the filter to see which one gives the highest psnr. The quality indices brisque, piqe and niqe were also used, but according to our analysis, only piqe is more coherent with the results and that is why its percentage of decrease is shown in ΔP (the smaller the piqe, the higher the image quality ).

- ![plot graphs image.](/Images/09.svg)

## GIFS

- **`get_gif_by_FWHM`**: inside the folder where the folders with the data are located, a folder will be created where inside there will be a gif with graphs and images that change according to the change in percent_FWHM, this type of gif is used to see the influence of this parameter In the images and curves associated with the experiment, all other parameters are maintained and only it will change. The user can only use create_gif, or create_gif_with_ifft to see the inverse transform in the time of the peaks.

- ![plot graphs image.](/Images/10.gif)

 - **`get_gif_by_interp_method`**: a gif will be created to see the difference between the interpolation types: 'nearest', 'linear', 'makima', 'pchip', 'spline'. The user can only use create_gif, as well as use create_gif_with_ifft to see the inverse transform in the time of the peaks.

- ![plot graphs image.](/Images/11.gif)

 - **`get_gif_by_ratiotukey`**: a gif will be created with the variation of the ratio in 0.05 intervals of the Tukey window used to smooth the signal after the initial impulse cut. The user can only use create_gif, as well as use create_gif_with_ifft to see the inverse transform in the time of the peaks.

- ![plot graphs image.](/Images/12.gif)

 - **`get_gif_by_ratiowindow`**: a gif will be created with the variation of the ratio in intervals of 0.01 of the window chosen by the user responsible for processing the signal after interpolation. The user can only use create_gif, as well as use create_gif_with_ifft to see the inverse transform in the time of the peaks.

- ![plot graphs image.](/Images/13.gif)

 - **`get_gif_by_windows`**: a gif will be created with all the window possibilities provided by matlab, all of them will be used with a ratio of 1. This way the user can choose the best window that gives the best spectrum and/or the best image. The user can only use create_gif, as well as use create_gif_with_ifft to see the inverse transform in the time of the peaks.

- ![plot graphs image.](/Images/14.gif)


PS: The gifs made are at a high speed, to be able to pause and view the image more calmly, the user can open the gif with Windows Media Player and pause to analyze it in more detail.



