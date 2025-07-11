
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

def run_tgif(fitsfile, catalogfile, band='b3', regionlabel='w51e', fix_pos_idx=None, fitting_size_dict=None):

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
                                vmin=None, vmax=None, maximum_size=4, savefitsdir='tables/', label_fits=f'{regionlabel}_{band}',fix_pos_idx=fix_pos_idx, fitting_size_dict=fitting_size_dict)

if __name__ == "__main__":
    run_tgif(Path.w51e_b6_conv, Path.w51e_dendro_matched_catalog_new, band='b6', regionlabel='w51e_conv', fix_pos_idx=[2,19,46,85,86,87], fitting_size_dict={2:0.4, 19:0.4, 46:0.4, 85:0.4, 86:0.4, 87:0.4})
    run_tgif(Path.w51n_b6_conv, Path.w51n_dendro_matched_catalog_new, band='b6', regionlabel='w51n_conv', fix_pos_idx=[35,37,39], fitting_size_dict={35:0.4, 37:0.8, 39:0.8})
    #run_tgif(Path.w51e_b3_tt0, Path.w51e_dendro_matched_catalog_new, band='b3', regionlabel='w51e')
    #run_tgif(Path.w51e_b6_tt0, Path.w51e_dendro_matched_catalog_new, band='b6', regionlabel='w51e')
    #run_tgif(Path.w51n_b3_tt0, Path.w51n_dendro_matched_catalog_new, band='b3', regionlabel='w51n')
    #run_tgif(Path.w51n_b6_tt0, Path.w51n_dendro_matched_catalog_new, band='b6', regionlabel='w51n')