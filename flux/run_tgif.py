
import sys
sys.path.append('/red/adamginsburg/t.yoo/anaconda3/lib/python3.9/site-packages/')
sys.path.append('/home/t.yoo/w51/TGIF')
#print(sys.path)

import TGIF.TGIF as tgif
import Paths.Paths as paths
from radio_beam import Beam
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table

import numpy as np
Path = paths.filepaths()

def run_tgif(fitsfile, catalogfile, band='b3', regionlabel='w51e'):

    fitsdata = fits.open(fitsfile)
    image = fitsdata[0].data
    if len(image.shape)>2:
        image= fitsdata[0].data[0][0]
   
    hdr = fits.getheader(fitsfile)  
    wcs = WCS(hdr,naxis=2)

    beam = Beam.from_fits_header(hdr)

    pixel_scale = wcs.proj_plane_pixel_scales()[0]

    catalog = Table.read(catalogfile)
    peakxy = np.vstack((catalog[f'{band}_xpix'], catalog[f'{band}_ypix'])).T

    peakxy_sky = np.vstack((catalog['b3_xsky'], catalog['b3_ysky'])).T


    tgif.plot_and_save_fitting_results(image, peakxy, beam, wcs, pixel_scale, fitting_size_default=0.6, saveimgdir=f'pngs/{regionlabel}_{band}/',label_img=f'{regionlabel}_{band}',
                                   vmin=None, vmax=None, maximum_size=4, savefitsdir='tables/', label_fits=f'{regionlabel}_{band}',)

if __name__ == "__main__":
    run_tgif(Path.w51e_b3_tt0, Path.w51e_dendro_matched_catalog_new, band='b3', regionlabel='w51e')
    run_tgif(Path.w51e_b6_tt0, Path.w51e_dendro_matched_catalog_new, band='b6', regionlabel='w51e')
    run_tgif(Path.w51n_b3_tt0, Path.w51n_dendro_matched_catalog_new, band='b3', regionlabel='w51n')
    run_tgif(Path.w51n_b6_tt0, Path.w51n_dendro_matched_catalog_new, band='b6', regionlabel='w51n')