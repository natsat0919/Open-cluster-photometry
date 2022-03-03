import os
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
from .a345_utilities import print_header, mem_usage
from matplotlib.patches import Rectangle as rect
from matplotlib.colors import LogNorm

class Correcting:
    def __init__(self, img_size_dir : str, target : str, raw_img_name : str):
        self.img_size_dir = img_size_dir
        self.target = target
        self.raw_img_name = raw_img_name

    def correcting(self):

        data_directory = f'/data/observatory/remote_telescope/data/cropped_{self.img_size_dir}'
        
        calibration_directory = os.path.join(data_directory, 'calibration/neg5c/master')  # neg5c refers to the CCD temperature during the calibration run 
        image_directory       = os.path.join(data_directory, f'images/Cluster Photometry/{self.target}')  # location of telescope image to be corrected
        
        image_filename = os.path.join(image_directory,f'{self.raw_img_name}')
        
        #load image to correct:
        with fits.open(image_filename) as hdu:
            image_header = hdu[0].header
            image_data = hdu[0].data
        
        exp_time = int(image_header["EXPTIME"])
        filtr = image_header["FILTER"][0]
        
        bias_filename = os.path.join(calibration_directory,'bias_master.fits')
        dark_filename = os.path.join(calibration_directory,f'dark_{exp_time}s_master.fits')
        flat_filename = os.path.join(calibration_directory,f'flat_{filtr}_master.fits')   ###XXX need proper image here

        

        # load each image and header
        with fits.open(bias_filename) as hdu:
            bias_header = hdu[0].header
            bias_data = hdu[0].data
        with fits.open(dark_filename) as hdu:
            dark_header = hdu[0].header
            dark_data = hdu[0].data
        with fits.open(flat_filename) as hdu:
            flat_header = hdu[0].header
            flat_data = hdu[0].data


        image_data_corrected = (image_data - dark_data)   / (flat_data - bias_data) * np.mean(flat_data - bias_data)


        #look at a small area of the image
        #image_data_corrected_sub = image_data_corrected[4000:5000, 4000:5000]
        #image_data_sub = image_data[4000:5000, 4000:5000]
        d_mean = np.mean(image_data);     # mean intensity
        d_std  = np.std(image_data);      # standard deviation of intensity
        d_min = np.min(image_data);
        d_max = np.max(image_data);
        d_med = np.median(image_data);
        d_lower=d_min
        d_upper=(d_max*d_med)**0.5
        
        d_mean_corrected = np.mean(image_data_corrected);     # mean intensity
        d_std_corrected  = np.std(image_data_corrected);      # standard deviation of intensity
        d_min_corrected = np.min(image_data_corrected);
        d_max_corrected = np.max(image_data_corrected);
        d_med_corrected = np.median(image_data_corrected);
        d_lower_corrected=d_med_corrected                     #good if we want to ignore background
        d_upper_corrected=(d_max_corrected*d_med_corrected)**0.5
        
        # plot the full corrected image
        #fig1, ax1 = plt.subplots(figsize=(18,20))
        #ax1.imshow(image_data, norm=LogNorm(vmin = d_lower, vmax=d_upper), cmap='gray')
        #ax1.set_title('Raw - full image')
        
        #fig2, ax2 = plt.subplots(figsize=(18,20))
        #ax2.imshow(image_data_corrected, norm=LogNorm(vmin = d_lower_corrected, vmax=d_upper_corrected), cmap='gray')
        #ax2.set_title('Corrected - full image')
        
        corrected_file_name = f"Corrected_{self.img_size_dir}_{self.raw_img_name[:-5]}"
        #plt.savefig(f"/data/observatory/student_data/natalia_bajnokova/{self.target}/Calibrated_{self.img_size_dir}/{filtr}/{corrected_file_name}.png")
        
        
        #Save to fits
        #hdu = fits.PrimaryHDU(data = image_data_corrected, header= image_header)
        #hdu.writeto(f"/data/observatory/student_data/natalia_bajnokova/Calibrated_{self.img_size_dir}/{self.filtr}/{corrected_file_name}.fits")
        
        
        return {"header":image_header, "data":image_data_corrected, "file_name":corrected_file_name, "img_size":self.img_size_dir, "target":self.target}
    
        #Plotting histogram
        #plt.figure(figsize=(9,6))
        #plt.hist(image_data_corrected.flatten(), bins=1000, range=(-2500, 2000), fc='k', ec='k');
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.xlabel('photo-electron counts in a pixel')
        #plt.ylabel('frequency')

        

    def calibration_imgs_analyse(self):
        data_directory = f'/data/observatory/remote_telescope/data/{self.img_size_dir}'
        calibration_directory = os.path.join(data_directory, 'calibration/neg5c/master')  # neg5c refers to the CCD temperature during the calibration run 

        bias_filename = os.path.join(calibration_directory,'bias_master.fits')
        dark_filename = os.path.join(calibration_directory,f'dark_{self.exp_time}s_master.fits')
        flat_filename = os.path.join(calibration_directory,f'flat_{self.filtr}_master.fits')   ###XXX need proper image here

        # load each image and header
        with fits.open(bias_filename) as hdu:
            bias_header = hdu[0].header
            bias_data = hdu[0].data
        with fits.open(dark_filename) as hdu:
            dark_header = hdu[0].header
            dark_data = hdu[0].data
        with fits.open(flat_filename) as hdu:
            flat_header = hdu[0].header
            flat_data = hdu[0].data

        #Mean and Standard deviations
            
        bias_mean = np.mean(bias_data)
        dark_mean = np.mean(dark_data)
        flat_mean = np.mean(flat_data)
        
        bias_std = np.std(bias_data)
        dark_std = np.std(dark_data)
        flat_std = np.std(flat_data)
    
        # Plot the calibration images, limiting the colour range to show detail. Also plot a small area zoomed in to show the pixel level variations    
        bias_min = bias_mean - 3*bias_std; bias_max = bias_mean + 3*bias_std
        dark_min = dark_mean - 3*dark_std; dark_max = dark_mean + 3*dark_std
        flat_min = flat_mean - 3*flat_std; flat_max = flat_mean + 3*flat_std
        nbins = 40

        f1, (ax_bias, ax_dark, ax_flat) = plt.subplots(1, 3, sharey=True, figsize=(20,8))
        plot_bias = ax_bias.imshow(bias_data, vmin = bias_min, vmax=bias_max, cmap='gray')
        plot_dark = ax_dark.imshow(dark_data, vmin = dark_min, vmax=dark_max, cmap='gray')
        plot_flat = ax_flat.imshow(flat_data, vmin = flat_min, vmax=flat_max, cmap='gray')
        ax_bias.set_title('Bias frame')
        ax_dark.set_title('Dark frame')
        ax_flat.set_title('Flat frame')

        #p=[1500,2000]                   
        #w=100; h=60;
        #ax_bias.add_patch(rect(p,w,h, edgecolor = 'firebrick', fill=False))
        #ax_dark.add_patch(rect(p,w,h, edgecolor = 'firebrick', fill=False))
        #ax_flat.add_patch(rect(p,w,h, edgecolor = 'firebrick', fill=False))
        #f2, (ax_bias_sub, ax_dark_sub, ax_flat_sub) = plt.subplots(1, 3, sharey=True, figsize=(20,8))
        #ax_bias_sub.imshow(bias_data[p[1]:p[1]+h,p[0]:p[0]+w], vmin = bias_min, vmax=bias_max, cmap='gray', interpolation="none")
        #ax_dark_sub.imshow(dark_data[p[1]:p[1]+h,p[0]:p[0]+w], vmin = dark_min, vmax=dark_max, cmap='gray', interpolation="none")
        #ax_flat_sub.imshow(flat_data[p[1]:p[1]+h,p[0]:p[0]+w], vmin = flat_min, vmax=flat_max, cmap='gray', interpolation="none")

        f3, (ax_bias_hist, ax_dark_hist, ax_flat_hist) = plt.subplots(1, 3, sharey=True, figsize=(20,6))
        ax_bias_hist.hist(bias_data.flatten(), bins=nbins, range=(bias_min, bias_max), fc='k', ec='k');
        ax_dark_hist.hist(dark_data.flatten(), bins=nbins, range=(dark_min, dark_max), fc='k', ec='k');
        ax_flat_hist.hist(flat_data.flatten(), bins=nbins, range=(flat_min, flat_max), fc='k', ec='k');