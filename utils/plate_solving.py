#%matplotlib widget
from astropy.io import fits, ascii
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
import os, time
from photutils import IRAFStarFinder
from astropy.table import QTable
from astropy.stats import mad_std
from photutils import aperture_photometry, CircularAperture
from .a345_utilities import print_header         # to use this in your own scripts copy the file a345_utilities.py to your work directory
import astropy.units as unit
from astropy.wcs.utils import pixel_to_skycoord

#set the default plot size 
plt.rcParams['figure.figsize'] = [15,15] # set the image size


# plate solver settings
print('Plate Solving module - Ignore any API key warnings below, they are not important')
from astroquery.astrometry_net import AstrometryNet
from astroquery.exceptions import TimeoutError
ast = AstrometryNet()
ast.API_URL = 'http://nova.astro.gla.ac.uk/api' # local server
ast.api_key = 'XXXXXXXX'
ast.URL = 'http://nova.astro.gla.ac.uk'

# Backup server on Ettus - use this as a back-up if there are issues with the primary server 
#ast.API_URL = 'http://ettus3.astro.gla.ac.uk:8000/api' # local server
#ast.URL = 'http://ettus3.astro.gla.ac.uk:8000'
#ast.api_key = 'XXXXXXXX'

print('')

def platesolving(**kwargs):

  
    #Load in the image
    #with fits.open(filename) as hdu:
    #    header = hdu[0].header
    #    data = hdu[0].data
    """
    header = image_details[0]
    data = image_details[1]
    corrected_img_name = image_details[2]
    img_size_dir = image_details[3]
    target = image_details[4]
    filtr = header["FILTER"][0] 
    """
    header = kwargs["header"]
    data = kwargs["data"]
    corrected_img_name = kwargs["file_name"]
    img_size_dir = kwargs["img_size"]
    target = kwargs["target"]
    filtr = header["FILTER"][0] 
    
    #some arrays
    unknown_orintation = []
    
    # get some colour scaling for the image
    d_mean = np.mean(data);     # mean intensity
    d_std  = np.std(data);      # standard deviation of intensity
    vmin = d_mean - d_std/2     # brightness at colour scale minimum
    vmax = d_mean + d_std       # brightness at colour scale maximum
    
    #rotate data
    if corrected_img_name[-6].lower() == "w":  
        data = np.rot90(data,1)
        y = True
    elif corrected_img_name[-6].lower() == "e":
        data = np.rot90(data,3)
        y = False
        
    else:
        unknown_orintation.append(corrected_img_name)
        az = header['OBJCTAZ'].split(' ')
        if 0 <= int(az[0]) <= 180:
            data = np.rot90(data,1)
            y= False
        else:
            data = np.rot90(data,3)
            y = True

    #Table of sources
    
    bkg_sigma = mad_std(data)                                    # get a measure of the image noise level
    NovaStars = IRAFStarFinder(fwhm=6, 
                               threshold=20.*bkg_sigma, 
                               minsep_fwhm=3, 
                               brightest=100)
    
    NovaSources = NovaStars(data)  
    for col in NovaSources.colnames:  
        NovaSources[col].info.format = '%.6g'  # for consistent table output 
    
    
    #Plate Solving
    
    image_width, image_height = data.shape
    wcs_header = None             # this variable will hld the result from the solver when it has completed
    t_start = time.time()
    try:
        print('Sending data to AstrometryNet server:')

        wcs_header = ast.solve_from_source_list(NovaSources['xcentroid'], NovaSources['ycentroid'],
                                                image_width, image_height,
                                                solve_timeout=300)
        if wcs_header:   # This will be true (ie not 'None') if the image solves sucessfully, and will contain a fits header with the solved image parameters (WCS)
            print('\n -> Success. Solving took {:0.1f}s'.format(time.time()-t_start))
            #print_header(wcs_header)         # the print_header function is in a345_utilities - see import statement at start. You can just use print(header),
        else:
            print('\n -> Solving failed after {:0.1f}s'.format(time.time()-t_start))
    except TimeoutError:
        #If you run out of time         
        print('\n -> ##FAIL: Timeout while solving, try a longer timeout, optmise number of sources (200-800 seems about right)')
        
        
    # update the existing header with the WCS data:
    print('Updated header, with WCS added')
    header.update(wcs_header)   # add in the solved WCS info
    
    wcs = WCS(header)
    
    # The easiest way to save the fils file is to create a new  fits abject, and add the data and calibrated header 
    hdu = fits.PrimaryHDU()         # create a FITS HDU object
    hdu.header.update(header)       # add in the header containg the wcs data
    hdu.data = data                 # add in the image data 
    
    # plot the fit image to save as png
    png_file = f'/data/observatory/student_data/natalia_bajnokova/{target}/Calibrated_{img_size_dir}/{filtr}/PlateSolved_{corrected_img_name}.png'
    ax = plt.subplot(projection=wcs)
    
    #get some colour limits
    rng=1.2
    cmin = d_mean/rng; cmax = d_mean*rng;
    
    # flip the y-axis to get dec increasing to the top
    if y:
        plt.gca().invert_yaxis()
    else:
        plt.gca().invert_xaxis()
    plt.imshow(data,vmin=cmin, vmax = cmax, cmap='gray')    
    plt.grid(color='lightblue', ls='solid')
    
    plt.savefig(png_file)
    plt.close()

    output_file = f'/data/observatory/student_data/natalia_bajnokova/{target}/Calibrated_{img_size_dir}/{filtr}/PlateSolved_{corrected_img_name}.fits'
    print('Saving calibrated file: {}'.format(output_file))
    hdu.writeto(output_file, overwrite=True)
    
    

    """
    Platesolved stars table
    """
    #Create table of RA, dec, flux, mag
    flux = NovaSources["flux"]
    peak = NovaSources["peak"]
    
    #Convert pixels to RA and dec
    i, sky = zip(*enumerate(wcs.pixel_to_world(NovaSources["xcentroid"],NovaSources["ycentroid"])))
    ra = [s.ra for s in sky]
    dec = [s.dec for s in sky]

    #Assign converted coordinates
    #ra = [c.ra.to_string(unit = unit.hour, sep=("h ", "m ", "s ")) for c in coord]
    #dec = [c.dec.to_string(unit = unit.degree, sep=("d ", "m ", "s ")) for c in coord]

    #create table
    
    t = QTable([i, ra , dec , flux, peak], names=("Index", 'RA','DEC', 'flux', "peak"), masked=True)
    
    for row in t:
        if row["peak"] > 60000:
            row["flux"] = "NaN"
        pass
    t.sort("flux")
    ascii.write(t, f'/data/observatory/student_data/natalia_bajnokova/{target}/Calibrated_{img_size_dir}/{filtr}/data_{corrected_img_name}.dat')
    
    return header, data