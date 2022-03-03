from astropy.io import fits, ascii
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
import os, time
from photutils import DAOStarFinder
from astropy.table import QTable
from astropy.stats import mad_std
from photutils import aperture_photometry, CircularAperture, CircularAnnulus
import astropy.units as unit
from astropy.wcs.utils import pixel_to_skycoord
from photutils import CircularAperture, CircularAnnulus, SkyCircularAperture, aperture_photometry
from matplotlib.colors import LogNorm
from photutils import IRAFStarFinder
from astropy.stats import sigma_clipped_stats
from os import listdir


class Photometry:
    
    @staticmethod
    def get_photometry(src, size, folderpath, stars_table):

        ra = stars_table["RA"]
        dec = stars_table["DEC"]
        filenames = [f for f in listdir(folderpath) if f.endswith('.fits')]
        filenames.sort()

        stars = {}
        for file in (filenames):

            print (f"Processing photometry on {file}")

            #import fit
            with fits.open(folderpath+file) as hdu:
                header = hdu[0].header
                data = hdu[0].data

            filtr = header["FILTER"][0]
            wcs = WCS(header)

            bkg_sigma = mad_std(data)        # get a measure of the image noise level

            SeeingStars = IRAFStarFinder(fwhm=8, threshold=50*bkg_sigma, brightest=250)     # set the detection threshold for a source based on the image noise
            SeeingStars = SeeingStars(data)
            fwhm = np.mean(SeeingStars['fwhm']) # calculates the average size of the stars in the image

            # convert ra and dec positions to pixel positions in the image
            xc,yc = wcs.wcs_world2pix(ra, dec, 1)
            PixPos = np.transpose((xc,yc))

            # sets the aperature and radii to use
            pix_radius = 4*fwhm/2 # 4 times the average star raduis to collect all the possible photons
            In_ann = pix_radius*1.2; Out_ann = (In_ann+1)*1.3 # annulus measurements from 1.2*aperture to 1.3*aperture

            # aperture photometry data
            aperture = CircularAperture(PixPos, r=pix_radius) # creating the cicular apertures cnetred on the star positions
            annulus_aperture = CircularAnnulus(PixPos, r_in=In_ann, r_out=Out_ann) # creating the annulus apertures cnetred on the star positions
            apers = [aperture, annulus_aperture]

            # Send the created apertures and the calibrated data to do the photometry
            print('stars and apertures created with wcs_world2pix, starting aperture photometry')
            phot_table = aperture_photometry(data, apers)

            phot_table['RA'] = ra; phot_table['DEC'] = dec                  # add to the table columns for Ra and Dec coordinates of the stars
            bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area # average value of the annulus around each star
            phot_table['background_mean'] = bkg_mean                        # adding the background level to the output table
            bkg_sum = bkg_mean * aperture.area                              # background for each star
            final_sum = phot_table['aperture_sum_0'] - bkg_sum              # The actual photometry results
            filter_counts_name = 'counts_'+header['FILTER']
            filter_flux_name = 'flux_'+header['FILTER']
            phot_table[filter_counts_name] = final_sum                      #  adding results to the output table
            filter_flux = final_sum/header['EXPTIME']
            phot_table[filter_flux_name] = filter_flux
            phot_table[f"mag {header['FILTER']}"] =-2.5*np.log10(filter_flux)
            
            ra_stars, dec_stars = wcs.wcs_pix2world(phot_table["xcenter"],phot_table["ycenter"], 1)
            
           
            flux = phot_table[filter_flux_name]
            count = phot_table[filter_counts_name]
            mag = phot_table[f"mag {header['FILTER']}"]
            tabl = QTable([ra_stars , dec_stars, flux, count, mag], names=('RA','DEC', "flux", "count", "mag"), masked=True)
            tabl.sort("flux")
            tabl.reverse()

            ascii.write(tabl, f'/data/observatory/student_data/natalia_bajnokova/{src}/Calibrated_{size}/{filtr}/data_{file}.dat')
            
            #stars = Photometry.stars_list(tabl, stars, filtr)
            
        #return stars
            
            
    @staticmethod
    def stars_list(tbl, stars, filtr):
        for row in tbl:
            if row["mag"] == "nan":
                pass
            else:
                r = round(row["RA"], 4)
                d = round(row["DEC"], 4)
                k = f"{r}, {d}"
                if k in stars:
                    if filtr == "G":
                        obj = stars[k]
                        obj.g.add(row["flux"], row["mag"])
                    elif filtr == "I":
                        obj = stars[k]
                        obj.i.add(row["flux"], row["mag"])
                    else:
                        obj = stars[k]
                        obj.r.add(row["flux"], row["mag"])
                else:
                    if filtr == "G":
                        stars[k] = Star(r, d)
                        stars[k].g.add(row["flux"], row["mag"])
                    elif filtr == "I":
                        stars[k] = Star(r, d)
                        stars[k].i.add(row["flux"], row["mag"])
                    else:
                        stars[k] = Star(r, d)
                        stars[k].r.add(row["flux"], row["mag"])
        return stars
        
    @staticmethod
    def average_magnitudes():
        pass
        
    
class Star:
    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec
        self.r = Filter()
        self.g = Filter()
        self.i = Filter()

    

class Filter:
        
    def __init__(self): 
        self.flux = []
        self.mag = []
        self.avg_mag = 0
        self.avg_flux = 0
            
    def add(self, flux, mag):
        self.flux.append(flux)
        self.mag.append(mag)
        
        if not np.isnan(mag) and not np.isnan(flux):
            self.avg_mag = np.mean(self.mag)
            self.avg_flux = np.mean(self.flux)